from time import time

def timeme(func):
    # This function shows the execution time of 
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = time()
        result = func(*args, **kwargs)
        t2 = time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result
    return wrap_func

def format_bytes(size: int, base: int = 2) -> str:
    if base not in (2,10):
        raise ValueError(f'Can only deal with base 2 or 10')
    # 2**10 = 1024
    power = 2**10 if base == 2 else 10**3
    n = 0
    power_labels = { 2: {0 : '', 1: 'Ki', 2: 'Mi', 3: 'Gi', 4: 'Ti'},
                    10: {0 : '', 1: 'K', 2: 'M', 3: 'G', 4: 'T'}}
    while size > power:
        size /= power
        n += 1
    return f'parsed {size} {power_labels.get(base)[n]}B'