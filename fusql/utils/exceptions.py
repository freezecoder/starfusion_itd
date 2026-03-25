"""Custom exceptions for FusionSQL."""


class FusionSQLError(Exception):
    """Base exception for FusionSQL errors."""
    pass


class ParserError(FusionSQLError):
    """Exception raised for parsing errors."""
    pass


class LoaderError(FusionSQLError):
    """Exception raised for data loading errors."""
    pass


class FileDiscoveryError(FusionSQLError):
    """Exception raised for file discovery errors."""
    pass
