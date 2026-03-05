class MoleculeError(Exception):
    """Base class for all errors in this package."""
    pass


class FragmentNotFoundError(MoleculeError):
    """Raised when a fragment is not found."""
    pass


class ManyFragmentsFoundError(MoleculeError):
    """Raised when too many fragments are found."""
    pass


class ChargeError(MoleculeError):
    """Raised when a charge of a molecule is invalid."""
    pass


class MultiplicityError(MoleculeError):
    """Raised when a spin multiplicity of a molecule is invalid."""
    pass