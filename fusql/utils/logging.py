import sys
from loguru import logger
from typing import Optional


def setup_logging(level: str = "INFO") -> None:
    """Configure loguru logging with the specified level.
    
    Args:
        level: Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    """
    logger.remove()
    logger.add(
        sys.stderr,
        level=level,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
        colorize=True,
    )


def get_logger(name: Optional[str] = None):
    """Get a logger instance.
    
    Args:
        name: Optional logger name for module identification
        
    Returns:
        Logger instance from loguru
    """
    if name:
        return logger.bind(name=name)
    return logger
