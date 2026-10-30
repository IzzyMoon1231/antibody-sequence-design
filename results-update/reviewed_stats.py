import statistics

def stats(xs):
 return {'n':len(xs),'mean':statistics.fmean(xs),'median':statistics.median(xs),'min':min(xs),'max':max(xs)}
