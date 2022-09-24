
def setup_args_log(ap, **args):
    import logging
    name = "Logger"
    if args.get("name"):
        name = args["name"]
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    return ap

def setup_logging(**args):
    """TODO: return logger as args["logger"]"""
    return args
