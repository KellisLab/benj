#import logging
class stopwatch:
    """usage:
        swgen = stopwatch.template("[INTEGRATION]")

        ...

        with swgen("Running xxx") as _:
            run_stuff()
        with swgen("Finalizing xxx") as _:
            finish_stuff()
        """
    def __init__(self, message, logger):
        self.logger = logger
        self.pre_message = message
        if len(message) > 1:
            self.post_message = message[0].lower() + message[1:]
        else:
            self.post_message = message
    def __enter__(self):
        from time import time
        self.logger.info(self.pre_message)
        self.timer = time()
        return self
    def tqdm_range(self, item_list, **kwargs):
        from tqdm.auto import tqdm
        return tqdm(item_list, desc=self.pre_message, **kwargs)
    def tqdm(self, **kwargs):
        return tqdm.tqdm(desc=self.pre_message, **kwargs)
    def __exit__(self, exc_type, exc_val, exc_tb):
        from time import time
        delta = time() - self.timer
        self.logger.info("Finished %s in %.2f seconds" % (self.post_message, delta))

def template(logname : str, level=None):
    import logging
    logger = logging.getLogger(logname)
    if level is not None:
        logger.setLevel(level)
    else:
        logger.setLevel(logging.INFO)
    return lambda msg: stopwatch(msg, logger=logger)
