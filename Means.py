## Written by Maxie D. Schmidt on 2017.06.19

from sage.all import *
from operator import mul

class MeanStats(object): 

     @staticmethod
     def power_mean(numlst, ppow): 
          nsize = float(len(numlst))
          eltpowsum = sum(map(lambda xi: (xi**ppow) / nsize, numlst))
          return eltpowsum.nth_root(ppow)
     ##
     
     @staticmethod
     def carleman1(numlst): ## See: http://mathworld.wolfram.com/CarlemansInequality.html
          product_func = lambda lst: reduce(mul, lst)
          prodlst = []
          for i in range(0, len(numlst)): 
               ithprod = product_func(numlst[0:i+1])
               prodlst += [ithprod.nth_root(i)]
          ## 
          return sum(prodlst)
     ## 
     
     @staticmethod
     def carleman2(numlst, ppow = 2.0): 
          sumlst = []
          for i in range(0, len(numlst)): 
               ithterm = float(sum(numlst[0:i+1]) / i)**ppow
               sumlst += [ithterm]
          ##
          return sum(sumlst)
     ## 
     
     @staticmethod
     def manhattan(numlst): 
          return sum(map(lambda xi: abs(xi), numlst))
     ## 
     
     @staticmethod
     def absolute_deviation(numlst): 
          mean = scipy.stats.mean(numlst)
          return sum(map(lambda xi: abs(xi - mean), numlst))
     ##
     
     @staticmethod
     def mean_deviation(numlst): 
          return float(Means.absolute_deviation(numlst) / len(numlst))
     ##
     
     @staticmethod
     def pearson_mode_skewness(numlst): ## See: http://mathworld.wolfram.com/PearsonModeSkewness.html
          mean = scipi.stats.mean(numlst)
          mode = scipy.stats.mode(numlst)
          stddev = scipi.stats.tstd(numlst)
          return float((mean - mode) / stddev)
     ## 
          
     @staticmethod
     def cubic_mean(numlst): 
          return Means.power_mean(numlst, ppow = 3)
     ##
     
     @staticmethod
     def contraharmonic(numlst): ## See: https://en.wikipedia.org/wiki/Contraharmonic_mean
          cnum = Means.power_mean(numlst, ppow = 2)**2
          cdenom = Means.power_mean(numlst, ppow = 1)
          return float(cnum / cdenom)
     ##
     
     @staticmethod
     def logphisum(numlst): 
          nsize = len(numlst)
          logphiterms = map(lambda xi: float(log(xi, golden_ratio) / nsize), numlst)
          return sum(logphiterms)
     ## 
     
     @staticmethod
     def coefficient_variation(numlst): ## See: https://en.wikipedia.org/wiki/Coefficient_of_variation
          mean = scipi.stats.mean(numlst)
          stddev = scipi.stats.tstd(numlst)
          return float(stddev / mean)
     ## 
     
     @staticmethod
     def interquartile_mean(numlst): ## See: https://en.wikipedia.org/wiki/Interquartile_mean
          nsize = len(numlst)
          sumlower, sumupper = floor(nsize / 4.0) + 1, floor(3.0 * nsize / 4.0)
          xiterms = numlst[sumlower:sumupper+1]
          return sum(map(lambda xi: 2.0 / nsize * xi, xiterms))
     ##
     
     @staticmethod
     def lmoment(numlst, moment = 2): ## See: https://en.wikipedia.org/wiki/L-moment
          nsize = len(numlst)
          lmomcoeff = lambda i: ([((-1.0)**k) * binomial(i-1, moment-1-k) * binomial(nsize-i, k) \
                                  for i in range(0, moment)])
          termlst = map(lambda (idx, xi): lmomcoeff(idx+1) * xi, enumerate(numlst))
          return float(sum(termlst) / binomial(nsize, moment) / moment)
     ##
     
     @staticmethod
     def compute_all_statistics(numlst): 
          means = [
               ("minimum", min), 
               ("maximum", max), 
               ("arithmetic", scipy.stats.mean), 
               ("geometric", scipy.stats.gmean), 
               ("harmonic", scipy.stats.hmean), 
               ("RMS", lambda lst: Means.power_mean(lst, ppow = 2)), 
               ("Carleman-v1", Means.carleman1), 
               ("Carleman-v2", Means.carleman2), 
               ("mode", scipy.stats.mode), 
               ("median", scipy.stats.median), 
               ("Manhattan", Means.manhattan), 
               ("variance", scipy.stats.tvar), 
               ("variation", scipy.stats.variation), 
               ("range", scipy.stats.range), 
               ("stddev", scipy.stats.tstd), 
               ("absdev", Means.absolute_deviation), 
               ("kurtosis", scipy.stats.mstats.kurtosis), 
               ("mean-deviation", Means.mean_deviation), 
               ("pearson-mode-skewness", Means.pearson_mode_skewness), 
               ("cubic", Means.cubic_mean), 
               ("interquartile-range", scipy.stats.iqr), 
               ("contraharmonic", Means.contraharmonic), 
               ("log-phi-sum", Means.logphisum), 
               ("coeff-of-variation", Means.coefficient_variation), 
               ("interquartile", Means.interquartile_mean), 
               ("lmoment2", lambda lst: Means.lmoment(lst, 2)), 
               ("lmoment3", lambda lst: Means.lmoment(lst, 3)), 
               ("lmoment4", lambda lst: Means.lmoment(lst, 4)), 
               ("circmean", scipy.stats.circmean), 
               ("circvar", scipy.stats.circmean), 
               ("circstd", scipy.stats.circstd), 
               ("skewness", scipy.stats.skew), 
               ("signal2noise", scipy.stats.signaltonoise)
          ] 
          rlst = [(mdesc, float(mfunc(numlst))) for (desc, mfunc) in means]
          return rlst
     ##
               
               
