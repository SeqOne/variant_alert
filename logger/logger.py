import logging


def get_logger(scope: str, level=logging.DEBUG):
    """
    get_logger
    """
    logging.basicConfig(
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=level
    )
    return logging.getLogger(scope)
