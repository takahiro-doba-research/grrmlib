class MoleculeError(Exception):
    """Base class for all errors in this package."""
    pass


class FragmentNotFoundError(MoleculeError):
    """Raised when a fragment is not found."""
    pass