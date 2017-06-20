## Written by Maxie D. Schmidt on 2017.06.19

from sage.all import *
from operator import mul
import scipy
from scipy.stats import *
from scipy.stats.mstats import hdmedian
import numpy as np

class MeanStats(object): 

     @staticmethod
     def power_mean(numlst, ppow): 
          nsize = float(len(numlst))
          eltpowsum = sum(map(lambda xi: (xi**ppow) / nsize, numlst))
          return eltpowsum**(1.0 / ppow)
     ##
     
     @staticmethod
     def carleman1(numlst): ## See: http://mathworld.wolfram.com/CarlemansInequality.html
          product_func = lambda lst: reduce(mul, lst)
          prodlst = []
          for i in range(1, len(numlst) + 1): 
               ithprod = product_func(numlst[0:i])
               prodlst += [ithprod**(1.0 / i)]
          ## 
          return sum(prodlst)
     ## 
     
     @staticmethod
     def carleman2(numlst, ppow = 2.0): 
          sumlst = []
          for i in range(1, len(numlst) + 1): 
               ithterm = float(sum(numlst[0:i]) / i)**ppow
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
          mean = scipy.stats.tmean(numlst)
          return sum(map(lambda xi: abs(xi - mean), numlst))
     ##
     
     @staticmethod
     def mean_deviation(numlst): 
          return float(MeanStats.absolute_deviation(numlst) / len(numlst))
     ##
     
     @staticmethod
     def pearson_mode_skewness(numlst): ## See: http://mathworld.wolfram.com/PearsonModeSkewness.html
          mean = scipy.stats.tmean(numlst)
          mode = float(scipy.stats.mode(numlst, axis = None).mode[0])
          stddev = 1 if scipy.stats.tstd(numlst) == np.nan else scipy.stats.tstd(numlst)
          return float((mean - mode) / stddev)
     ## 
          
     @staticmethod
     def cubic_mean(numlst): 
          return MeanStats.power_mean(numlst, ppow = 3)
     ##
     
     @staticmethod
     def contraharmonic(numlst): ## See: https://en.wikipedia.org/wiki/Contraharmonic_mean
          cnum = MeanStats.power_mean(numlst, ppow = 2)**2
          cdenom = MeanStats.power_mean(numlst, ppow = 1)
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
          mean = scipy.stats.tmean(numlst)
          stddev = 1 if scipy.stats.tstd(numlst) == np.nan else scipy.stats.tstd(numlst)
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
          if nsize < moment: 
               return 0.0
          ##
          lmomcoeff = lambda i: sum([((-1.0)**k) * binomial(i-1, moment-1-k) * binomial(nsize-i, k) \
                                  for k in range(0, moment)])
          termlst = map(lambda (idx, xi): lmomcoeff(idx+1) * xi, enumerate(numlst))
          return float(sum(termlst) / binomial(nsize, moment) / moment)
     ##
     
     @staticmethod
     def compute_all_statistics(numlst): 
          means = [
               ("minimum", min), 
               ("maximum", max), 
               ("arithmetic", scipy.stats.tmean), 
               ("geometric", scipy.stats.gmean), 
               ("harmonic", scipy.stats.hmean), 
               ("RMS", lambda lst: MeanStats.power_mean(lst, ppow = 2)), 
               #("Carleman-v1", MeanStats.carleman1), 
               #("Carleman-v2", MeanStats.carleman2), 
               ("mode", lambda lst: scipy.stats.mode(lst).mode[0]), 
               ("median", lambda lst: scipy.stats.mstats.hdmedian(lst).data), 
               #("Manhattan", MeanStats.manhattan), 
               #("variance", scipy.stats.tvar), 
               #("variation", scipy.stats.variation), 
               #("stddev", scipy.stats.tstd), 
               #("absdev", MeanStats.absolute_deviation), 
               #("kurtosis", scipy.stats.mstats.kurtosis), 
               #("mean-deviation", MeanStats.mean_deviation), 
               #("pearson-mode-skewness", MeanStats.pearson_mode_skewness), 
               #("cubic", MeanStats.cubic_mean), 
               #("interquartile-range", scipy.stats.iqr), 
               #("contraharmonic", MeanStats.contraharmonic), 
               #("log-phi-sum", MeanStats.logphisum), 
               #("coeff-of-variation", MeanStats.coefficient_variation), 
               #("interquartile", MeanStats.interquartile_mean), 
               #("lmoment2", lambda lst: MeanStats.lmoment(lst, 2)), 
               #("lmoment3", lambda lst: MeanStats.lmoment(lst, 3)), 
               #("lmoment4", lambda lst: MeanStats.lmoment(lst, 4)), 
               #("circmean", scipy.stats.circmean), 
               #("circvar", scipy.stats.circmean), 
               #("circstd", scipy.stats.circstd), 
               #("skewness", scipy.stats.skew), 
               ("signal2noise", lambda lst: float(scipy.stats.signaltonoise(lst)))
          ] 
          rlst = [(mdesc, mfunc(numlst)) for (mdesc, mfunc) in means]
          return rlst
     ##
               
               
