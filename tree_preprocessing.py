"""
Tree Preprocessing Module for TreeHarmonizer

This module handles the normalization of arbitrary newick trees to the format
required by TreeHarmonizer. The expected format includes:
1. An "extra root" node - a single-child root wrapping the entire tree
2. "Private nodes" - intermediate nodes between each leaf and its parent

Trees not in this format are automatically transformed.
"""

from ete3 import Tree


# Constants for naming synthetic nodes
ROOT_NAME = "__TH_ROOT__"
PRIVATE_SUFFIX = "_private"


class TreeMetadata:
    """
    Metadata about tree structure and naming mappings.

    This class tracks the relationship between original and transformed
    node names, enabling downstream code to work with arbitrary input trees.
    """

    def __init__(self, original_newick, was_transformed, leaf_names=None,
                 private_node_to_leaf=None, leaf_to_private_node=None,
                 root_name="", had_extra_root=False, had_private_nodes=False):
        """
        Initialize TreeMetadata.

        Args:
            original_newick: The original newick string before transformation
            was_transformed: Whether the tree was transformed
            leaf_names: List of leaf/sample node names
            private_node_to_leaf: Dict mapping private node names to leaf names
            leaf_to_private_node: Dict mapping leaf names to private node names
            root_name: Name of the root node
            had_extra_root: Whether the original tree had an extra root
            had_private_nodes: Whether the original tree had private nodes
        """
        self.original_newick = original_newick
        self.was_transformed = was_transformed
        self.leaf_names = leaf_names if leaf_names is not None else []
        self.private_node_to_leaf = private_node_to_leaf if private_node_to_leaf is not None else {}
        self.leaf_to_private_node = leaf_to_private_node if leaf_to_private_node is not None else {}
        self.root_name = root_name
        self.had_extra_root = had_extra_root
        self.had_private_nodes = had_private_nodes

    def is_leaf(self, node_name):
        """Check if node name is a leaf/sample node."""
        return node_name in self.leaf_names

    def is_private_node(self, node_name):
        """Check if node is a private node (directly above a leaf)."""
        return node_name in self.private_node_to_leaf

    def is_internal_node(self, node_name):
        """
        Check if node is a true internal node.

        Internal nodes are not: root, private nodes, or leaves.
        """
        return (
            not self.is_leaf(node_name) and
            not self.is_private_node(node_name) and
            node_name != self.root_name
        )

    def get_leaf_for_private(self, private_name):
        """Get the leaf name for a given private node."""
        return self.private_node_to_leaf.get(private_name)

    def get_private_for_leaf(self, leaf_name):
        """Get the private node name for a given leaf."""
        return self.leaf_to_private_node.get(leaf_name)


def has_extra_root(tree: Tree) -> bool:
    """
    Detect if tree has an "extra root" node.

    An extra root exists if:
    1. Root has exactly one child
    2. That child's descendant leaves equal the root's leaves

    Args:
        tree: An ete3 Tree object

    Returns:
        True if extra root exists, False if it needs to be added
    """
    root = tree.get_tree_root()
    children = root.get_children()

    # Must have exactly one child
    if len(children) != 1:
        return False

    # Child's leaves must equal root's leaves
    root_leaves = set(root.get_leaf_names())
    child_leaves = set(children[0].get_leaf_names())

    return root_leaves == child_leaves


def has_private_nodes(tree: Tree) -> bool:
    """
    Detect if tree has "private nodes" before each leaf.

    Private nodes exist if each leaf's parent:
    1. Has exactly one child (the leaf)
    2. Is not the root

    Args:
        tree: An ete3 Tree object

    Returns:
        True if all leaves have private nodes, False otherwise
    """
    root = tree.get_tree_root()
    leaves = tree.get_leaves()

    if not leaves:
        return False

    for leaf in leaves:
        parent = leaf.up

        # Leaf's parent must exist and not be root
        if parent is None or parent == root:
            return False

        # Parent must have exactly one child (the leaf)
        if len(parent.get_children()) != 1:
            return False

    return True


def detect_tree_structure(tree: Tree) -> dict:
    """
    Analyze tree structure and return a report.

    Args:
        tree: An ete3 Tree object

    Returns:
        Dictionary with structure analysis results
    """
    extra_root = has_extra_root(tree)
    private_nodes = has_private_nodes(tree)

    return {
        'has_extra_root': extra_root,
        'has_private_nodes': private_nodes,
        'leaf_names': tree.get_leaf_names(),
        'needs_transformation': not (extra_root and private_nodes)
    }


def add_extra_root(tree: Tree, root_name: str = ROOT_NAME) -> Tree:
    """
    Wrap tree with an extra single-child root node.

    Before: ((A,B)C)D;   (D is current root)
    After:  (((A,B)C)D)__TH_ROOT__;

    Args:
        tree: An ete3 Tree object (modified in place)
        root_name: Name for the new root node

    Returns:
        The modified tree
    """
    # Create new root
    new_root = Tree()
    new_root.name = root_name

    # Get current root and make it a child of new root
    current_root = tree.get_tree_root()

    # Detach current root and add as child of new root
    new_root.add_child(current_root)

    return new_root


def add_private_nodes(tree: Tree, suffix: str = PRIVATE_SUFFIX) -> tuple:
    """
    Insert private nodes between each leaf and its parent.

    Before: (A,B)C;
    After:  ((A)A_private,(B)B_private)C;

    Args:
        tree: An ete3 Tree object (modified in place)
        suffix: Suffix to add to leaf name for private node

    Returns:
        Tuple of (modified tree, mapping dict {private_name: leaf_name})
    """
    mapping = {}
    leaves = tree.get_leaves()

    root = tree.get_tree_root()

    for leaf in leaves:
        parent = leaf.up
        if parent is None:
            continue

        # Create private node name
        private_name = f"{leaf.name}{suffix}"
        mapping[private_name] = leaf.name

        # Check if parent already has only this leaf as child AND is not the root
        # (meaning private node structure already exists for this leaf)
        # If parent IS the root, we still need to add a private node
        if len(parent.get_children()) == 1 and parent != root:
            continue

        # Create new private node
        private_node = Tree()
        private_node.name = private_name

        # Detach leaf from parent
        leaf.detach()

        # Add leaf as child of private node
        private_node.add_child(leaf)

        # Add private node to original parent
        parent.add_child(private_node)

    return tree, mapping


def _build_metadata_from_tree(tree: Tree, original_newick: str,
                               was_transformed: bool,
                               had_extra_root: bool,
                               had_private_nodes: bool) -> TreeMetadata:
    """
    Build TreeMetadata from the current tree state.

    Args:
        tree: The (possibly transformed) tree
        original_newick: The original input newick string
        was_transformed: Whether any transformation was applied
        had_extra_root: Whether original tree had extra root
        had_private_nodes: Whether original tree had private nodes

    Returns:
        TreeMetadata instance
    """
    root = tree.get_tree_root()
    leaves = tree.get_leaves()

    # Build leaf list
    leaf_names = [leaf.name for leaf in leaves]

    # Build private node mappings
    private_node_to_leaf = {}
    leaf_to_private_node = {}

    for leaf in leaves:
        parent = leaf.up
        if parent is not None and parent != root:
            # Check if parent is a private node (has only this leaf as child)
            if len(parent.get_children()) == 1:
                private_node_to_leaf[parent.name] = leaf.name
                leaf_to_private_node[leaf.name] = parent.name

    return TreeMetadata(
        original_newick=original_newick,
        was_transformed=was_transformed,
        leaf_names=leaf_names,
        private_node_to_leaf=private_node_to_leaf,
        leaf_to_private_node=leaf_to_private_node,
        root_name=root.name,
        had_extra_root=had_extra_root,
        had_private_nodes=had_private_nodes
    )


def normalize_tree(newick_string: str, verbose: bool = False) -> tuple:
    """
    Main entry point: detect and transform tree as needed.

    This function analyzes the input tree and adds any missing structural
    elements (extra root, private nodes) required by TreeHarmonizer.

    Args:
        newick_string: Newick format tree string
        verbose: If True, print information about transformations

    Returns:
        Tuple of (normalized newick string, TreeMetadata object)
    """
    # Parse tree
    tree = Tree(newick_string, format=1)

    # Detect current structure
    structure = detect_tree_structure(tree)
    had_extra_root = structure['has_extra_root']
    had_private_nodes = structure['has_private_nodes']

    was_transformed = False

    if verbose:
        print(f"TreeHarmonizer: Analyzing input tree...")
        print(f"  - Detected {len(structure['leaf_names'])} leaf nodes (samples)")

    # Add private nodes first (before adding root, so tree structure is correct)
    if not had_private_nodes:
        tree, _ = add_private_nodes(tree)
        was_transformed = True
        if verbose:
            print(f"  - Added private nodes for variant placement")

    # Add extra root if needed
    if not had_extra_root:
        tree = add_extra_root(tree)
        was_transformed = True
        if verbose:
            print(f"  - Added root node: {ROOT_NAME}")

    if verbose:
        if was_transformed:
            print(f"  - Tree normalization complete.")
        else:
            print(f"  - Tree already in required format, no transformation needed.")

    # Build metadata
    metadata = _build_metadata_from_tree(
        tree,
        newick_string,
        was_transformed,
        had_extra_root,
        had_private_nodes
    )

    # Return normalized newick string and metadata
    normalized_newick = tree.write(format=1)

    return normalized_newick, metadata


def validate_tree_for_placement(newick_string: str) -> None:
    """
    Validate that a newick string is suitable for TreeHarmonizer.

    Raises ValueError with descriptive message if:
    - Tree cannot be parsed
    - Tree has no leaves
    - Leaf names contain reserved prefixes/suffixes

    Args:
        newick_string: Newick format tree string

    Raises:
        ValueError: If tree is invalid for placement
    """
    try:
        tree = Tree(newick_string, format=1)
    except Exception as e:
        raise ValueError(f"Failed to parse newick string: {e}")

    leaves = tree.get_leaves()

    if not leaves:
        raise ValueError("Tree has no leaf nodes (samples)")

    # Check for reserved naming patterns
    for leaf in leaves:
        if leaf.name.endswith(PRIVATE_SUFFIX):
            raise ValueError(
                f"Leaf name '{leaf.name}' ends with reserved suffix '{PRIVATE_SUFFIX}'. "
                "Please rename your samples to avoid this suffix."
            )
        if leaf.name == ROOT_NAME:
            raise ValueError(
                f"Leaf name '{leaf.name}' is a reserved name. "
                "Please rename this sample."
            )

    # Check all nodes have names
    for node in tree.traverse():
        if not node.name or node.name.strip() == "":
            raise ValueError(
                "Tree contains unnamed nodes. All nodes (root, internal, and leaves) "
                "must have names in the newick string."
            )
