from sage.all import *
from sage.plot.histogram import Histogram, histogram
from itertools import takewhile
from ulam_sets import compute_ulam_set
from sage.plot.colors import rainbow
import numpy as np
from sympy import fourier_transform, inverse_fourier_transform
from scipy.stats.mstats import mquantiles
from numpy import unique, arange

import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Tiling import *

def V(x, y): 
     return vector([x, y])
##

def max_norm(v): 
     return max(map(abs, list(v)))
##

def edist(v): # === Euclidean distance squared
     xycoords = list(v)
     comps = map(lambda xycoord: xycoord ** 2, xycoords)
     return sum(comps)
##

def compute_Sn_set(ulam_set, norm_func, alpha = 2.5714474995): 
     Sn = []
     for an in ulam_set: ## my guess at what the set S_N should look like in 2d: 
          sn = norm_func(alpha * an) - float(2*pi) * floor(norm_func(alpha * an / 2.0 / float(pi)))
          Sn += [sn]
     ## 
     return Sn
##

def save_ulam_set_image_v2(outfile, init_vectors, n = 100, a1a2 = [1, 1], 
                           norm_funcs = [(max_norm, "Max Norm"), (edist, "Euclidean Distance")]):  
           
     # plot options (change these as you will to explore): 
     thickness, sbase, aratio, ps = 2, float(golden_ratio), 'automatic', 15
     scale_options = [('linear', 10), ('semilogy', sbase), ('semilogx', sbase), ('loglog', sbase)]
     
     # compute the plots: 
     init_conds_len = len(init_vectors)
     image_graphics = []
     for (norm_func, nfunc_desc) in norm_funcs: 
          ulam_set = compute_ulam_set_v2(a1a2[0], a1a2[1], n, init_vectors, norm_func)
          print "ULAM SET: ", ulam_set, "\n"
          gplot = point(ulam_set, pointsize=ps, axes = False, axes_labels = None, gridlines = None)
          snset = compute_Sn_set(ulam_set, norm_func)
          g1, g2, g3, g4 = Graphics(), Graphics(), Graphics(), Graphics()
          snhist = g1 + histogram(snset, normed = True)
          snhist_normed = g2 + histogram(snset, normed = True)
          snhist_normed_bins100 = g3 + histogram(snset, bins = 100, normed = True)
          snhist_normed_cdf = g4 + histogram(snset, normed = True, cumulative = True)
          graphics_row = [gplot, snhist, snhist_normed, snhist_normed_bins100, snhist_normed_cdf]
          for sopt in []: #scale_options: 
               plot_title = "Ulam Set Plots: $%s$\n\nUlam Set Plot\nNorm = %s, Scale = %s" % ("[]", nfunc_desc, sopt)
               #gplot = plot(ulam_set, title = plot_title, scale = sopt[0], base = sopt[1])
               gplot = point(ulam_set, pointsize=ps, scale = sopt, axes = False, 
                             axes_labels = None, gridlines = None, aspect_ratio = 'automatic', 
                             title = plot_title)
               #gplot.save("test-%s.png" % nfunc_desc, axes = False, axes_labels = None)
          ##
          image_graphics += [graphics_row]
          
     ## 
     garray = graphics_array(image_graphics)
     plot_title_lst = map(lambda (x, y): r'$[^{% 3g}_{% 3g}]$' % (x, y), init_vectors)[0:init_conds_len]
     plot_title = "\{" + ', '.join(plot_title_lst) + ", \ldots\}"
     garray.show(title = plot_title, fontsize = 5, frame = True, typeset = 'latex', axes = False, axes_labels = ("", ""))
     garray.save(outfile, title = plot_title, fontsize = 14, axes = False, frame = True, gridlines = False, axes_labels = ("", ""), typeset = 'latex', figsize = [10, 10])
     
## 

def save_ulam_set_image(outfile, init_vectors, n = 100, a1a2 = [1, 1], 
                        norm_funcs = [(max_norm, "Max Norm")]):  
           
     # plot options (change these as you will to explore): 
     thickness, sbase, aratio, ps = 2, float(golden_ratio), 'automatic', 15
     scale_options = [('linear', 10), ('semilogy', sbase), ('semilogx', sbase), ('loglog', sbase)]
     
     # compute the plots: 
     init_conds_len = len(init_vectors)
     image_graphics = []
     ulam_set, ablattice = [], []
     for (norm_func, nfunc_desc) in norm_funcs: 
          #ulam_set = compute_ulam_set_v2(a1a2[0], a1a2[1], n, init_vectors, norm_func)
          ulam_set = compute_ulam_set(n, init_vectors)
          print "ULAM SET: [LEN = %d]" % len(ulam_set), sorted(list(ulam_set)), "\n"
          print "ULAM SET: ", ulam_set, "\n"
          nupper = 75 # n + 1 #N(log(n) / log(golden_ratio))
          ulamset2 = ulam_set.difference(set([(k*golden_ratio+1, golden_ratio+k) for k in range(1, nupper)]))
          ulamset2 = ulamset2.difference(set([(golden_ratio+k, k*golden_ratio+1) for k in range(1, nupper)]))
          ablattice = set([(a*golden_ratio+b, b*golden_ratio+a) for a in range(1, nupper, 2) for b in range(1, nupper, 2) \
                           if b*golden_ratio+a < -(a*golden_ratio+b)+(nupper+1)*(golden_ratio+1)])
          print "LATTICE DIFFERENCE: ", ablattice.difference(ulam_set), "\n"
          
          lucas = lambda k: Integer((5/2)*lucas_number1(k-1,1,-1)+(1/2)*lucas_number2(k-1,1,-1))
          lucasodd = [lucas(k) for k in range(0, nupper) if lucas(k) % 2 == 1] 
          fibodd = [fibonacci(k) for k in range(0, nupper) if fibonacci(k) % 2 == 1]
          abtuple_func = lambda a, b: (a*golden_ratio+b, b*golden_ratio+a)
          lcombos = [abtuple_func(l1, l2) for l1 in lucasodd for l2 in lucasodd]
          fcombos = [abtuple_func(f1, f2) for f1 in fibodd for f2 in fibodd]
          lfcombos = [abtuple_func(l1, f2) for l1 in lucasodd for f2 in fibodd]
          flcombos = [abtuple_func(f2, l1) for l1 in lucasodd for f2 in fibodd]
          tuple_combos = lcombos + fcombos + lfcombos + flcombos
          ulamset2 = ulamset2.difference(set(tuple_combos))
          
          lucasoddadded = [lucas(k) + lucas(m) for k in range(0, nupper) for m in range(0, nupper) if (lucas(k) + lucas(m)) % 2 == 1]
          fiboddadded = [fibonacci(k) + fibonacci(m) for k in range(0, nupper) for m in range(0, nupper) if (fibonacci(k) + fibonacci(m)) % 2 == 1]
          fladded = [fibonacci(k) + lucas(m) for k in range(0, nupper) for m in range(0, nupper) if (fibonacci(k) + lucas(m)) % 2 == 1]
          addedlst = sorted(list(unique(lucasoddadded + fiboddadded + fladded)))[0:n+1]
          #print addedlst, "\n", lucasoddadded, "\n", fiboddadded, "\n", fladded, "\n"
          
          afcombos = [abtuple_func(a1, f2) for a1 in addedlst for f2 in fibodd]
          facombos = [abtuple_func(f2, a1) for a1 in addedlst for f2 in fibodd]
          alcombos = [abtuple_func(a1, l2) for a1 in addedlst for l2 in lucasodd]
          lacombos = [abtuple_func(l2, a1) for a1 in addedlst for l2 in lucasodd]
          aacombos = [abtuple_func(a1, a2) for a1 in addedlst for a2 in addedlst]
          added_tuple_combos = afcombos + facombos + alcombos + lacombos + aacombos
          ulamset2 = ulamset2.difference(set(added_tuple_combos))
          print set(range(0, n**2)).difference(lucasodd + fibodd + lucasoddadded + fiboddadded + fladded), "\n"
          
          print "ULAM SET (INTERIOR): ", sorted(list(ulamset2)), "\n"
                   
          tuple_combos2 = set(tuple_combos + added_tuple_combos)
          tuple_combos2 = tuple_combos2.difference(ulam_set)
          print "UNUSED VECTORS (INTERIOR): ", sorted(list(tuple_combos2))[0:16], "\n"
          
          #return None

          gplot = point(ulam_set, pointsize=ps, axes = False, axes_labels = None, gridlines = None)
          gplot += line([(golden_ratio+m,m*golden_ratio+1) for m in range(1, nupper)])
          gplot += line([(n*golden_ratio+1,golden_ratio+n) for n in range(1, nupper)])
          redlines = [line([(n*golden_ratio+1,golden_ratio+n), (golden_ratio+n, n*golden_ratio+1)], rgbcolor = 'red') for n in range(1, nupper)]
          for rline in redlines: gplot += rline
          gplot += point([(a*golden_ratio+b, b*golden_ratio+a) for a in range(1, nupper, 2) for b in range(1, nupper, 2) if b*golden_ratio+a < -(a*golden_ratio+b)+(nupper+1)*(golden_ratio+1)], axes = False, axes_labels = None, gridlines = None, rgbcolor = 'limegreen', marker="s", ps = 45)
          image_graphics += [[gplot]]
          break
          gaphist = Graphics()
          gappoints = Tiling.compute_slope_gaps(list(ulam_set))
          #gappoints = map(lambda pt: pt / len(list(ulam_set)), gappoints)
          [minr, maxr] = mquantiles(np.array(gappoints), prob = [0.01, 0.9])
          gaphist += histogram(gappoints, normed = True, label = 'Slope Gap Distribution', range = [minr, maxr], bins = 500)
          
          agaphist = Graphics()
          agappoints = Tiling.compute_angle_gaps(list(ulam_set))
          #gappoints = map(lambda pt: pt / len(list(ulam_set)), gappoints)
          [minr, maxr] = mquantiles(np.array(agappoints), prob = [0.01, 0.9])
          agaphist += histogram(agappoints, normed = True, label = 'Angle Gap Distribution', range = [minr, maxr], bins = 500)
          
          pchist = Graphics()
          pcpoints = Tiling.compute_pc_edists(list(ulam_set))
          #pcpoints = map(lambda pt: pt / len(list(ulam_set)), pcpoints)
          [minr, maxr] = mquantiles(np.array(pcpoints), prob = [0.01, 0.9])
          pchist += histogram(pcpoints, normed = True, label = 'Slope Gap Distribution', range = [minr, maxr], bins = 500)
          image_graphics += [[gplot], [gaphist], [agaphist], [pchist]]     
     ## 
     garray = graphics_array(image_graphics)
     plot_title_lst = map(lambda (x, y): r'$[^{% 3g}_{% 3g}]$' % (x, y), init_vectors)[0:init_conds_len]
     plot_title = "\{" + ', '.join(plot_title_lst) + ", \ldots\}"
     garray.show(title = plot_title, fontsize = 5, frame = True, typeset = 'latex', axes = False, axes_labels = ("", ""))
     garray.save(outfile, title = plot_title, fontsize = 14, axes = False, frame = True, gridlines = False, axes_labels = ("", ""), typeset = 'latex', figsize = [10, 10])
     
     return ulam_set, ablattice
     
## 

def save_example_images(n = 10, a1a2 = [1, 1]):    
     
     (a, b) = map(float, a1a2)
     absuffix = "a%03db%03d" % (a, b)
     initial_vector_configs = [ ## examples from Jayadev's talk and in the article: 
          #[V(1, golden_ratio), V(0, 1)], 
          [V(1, golden_ratio), V(golden_ratio, 1)], 
          #[V(1, float(golden_ratio)), V(float(golden_ratio), 1)], 
          #[V(1, golden_ratio), V(1, 0)], 
          #[V(1, 0), V(0, 1)], 
          #[V(9, 0), V(0, 9), V(1, 13)], 
          #[V(2, 5), V(3, 1)], 
          #[V(1, 0), V(2, 0), V(0, 1)], 
          #[V(2, 0), V(0, 1), V(3, 1)], 
          #[V(1, 0), V(0, 1), V(2, 3)], 
          #[V(3, 0), V(0, 1), V(1, 1)], 
          #[V(1, 0), V(2, 0), V(0, 1)], 
          #[V(2, 0), V(3, 0), V(0, 1)], 
          #[V(1, 0), V(0, 1), V(6, 4)], 
          #[V(1, 0), V(0, 1), V(10, 9)], 
          #[V(1, 0), V(0, 1), V(10, 3)], 
          #[V(1, 3), V(3, 4)], 
          #[V(1, 0), V(1, 1)]
     ] 
     
     for (icidx, init_vectors) in enumerate(initial_vector_configs): 
          plot_suffix = "ulam-set" + "-N." + "%05d" % n + "-" + absuffix + "-v" + str(icidx + 1) + ".png"
          print "  => Saving image \"%s\" ... " % plot_suffix
          save_ulam_set_image(plot_suffix, map(tuple, init_vectors), n, a1a2)
     ## 
     
##

def generate_lincombo_comp_graphs(outfile_suffix, init_vectors, n, norm_func = max_norm): 
     
     garray_data = []
     for a1 in range(1, 5): 
          graphics_row = []
          for a2 in range(1, 5): 
               print "a1 / a2", [a1, a2]
               ulam_set = compute_ulam_set_v2(a1, a2, n, init_vectors, norm_func) 
               print "ULAM_SET: ", ulam_set
               plot_title = r"$a_1$ / $a_2$ = % 3g / % 3g" % (a1, a2)
               gplot = point(ulam_set, pointsize=2, axes = False, axes_labels = None, gridlines = None, 
                             title = plot_title)
               graphics_row += [gplot]
               #graphics_row += [ulam_set]
          ##
          garray_data += [graphics_row]
     ## 
     
     outfile = 'ulam-set-' + outfile_suffix + '.png'
     garray = graphics_array(garray_data)
     garray.show(fontsize = 5, frame = True, typeset = 'latex', axes = False, axes_labels = ("", ""))
     garray.save(outfile, fontsize = 14, axes = False, frame = True, gridlines = False, 
                 axes_labels = ("", ""), figsize = [10, 10])

##

def defint_func(m, n, alpha): 
     if m == 0 and n == 0: 
          return 0
     elif m == 0: 
          return sin(2 * alpha * n * pi) / 2.0 / alpha / n / pi
     elif n == 0: 
          return sin(2 * alpha * m * pi) / 2.0 / alpha / m / pi
     else: 
          return cos(alpha *(m+n) * pi) * sin(alpha * m * pi) * sin(alpha * n * pi) / (alpha**2) / m / n / (pi ** 2)
     ##
## 

def compute_2d_integral(n, init_vectors = [V(1, float(golden_ratio)), V(1, 0)]): 

     ulam_set = compute_ulam_set(n, map(tuple, init_vectors))
     alpha, x, y = var('alpha x y')
     defint = sum(map(lambda (m, n): defint_func(float(m), float(n), alpha), ulam_set))
     return defint
     
     #eqns = [cos(alpha * (m+n)) < 0 for (m, n) in ulam_set]
     #solve(eqns, alpha)
     
     N = len(ulam_set)
     integrand = lambda x, y: sum(map(lambda (m, n): cos(2 * pi * alpha *(m*x+n*y)), ulam_set))
     defintx = lambda y: integral(integrand(x, y), x, 0, 1)
     defint = integral(defintx(y), y, 0, 1)
     graphics = plot(defint, (-2*pi, 2*pi))
     graphics += plot(N, (-2*pi, 2*pi), title = "N = %d" % N)
     plot_func = lambda alpha, beta: defint(alpha)-beta*N
     return defint, plot_func, graphics
     #print simplify(defint)
     #print find_root(defint == -0.8 * N, 0, 2*pi)
     #print defint.subs(alpha == 1.0).n()
     #print map(lambda soln: soln.rhs(), solve([defint == -0.8 * N], alpha))

##

def compute_2d_integral_plots_old(init_vectors = [V(1, golden_ratio), V(1, 0)]): 
     nvalues = [5, 10, 15, 20, 25, 30, 35, 40, 45, 45, 50, 55]
     garray, grow = [], []
     for (nidx, N) in enumerate(nvalues): 
          grow += [compute_2d_integral(N, init_vectors)]
          if nidx % 3 == 2: 
               garray += [grow]
               grow = []
          ##
     ##
     gplots = graphics_array(garray)
     gplots.show(frame = True)
##

def compute_2d_integral_plots(init_vectors = [V(1, golden_ratio), V(1, 0)], 
                              plot_id = "V2"): 
     nvalues = [25, 50, 100, 150, 250, 350, 500, 650, 750, 850, 950, 1050, 1250, 1500, 1600, 1700, 1800, 1900, 2000, 2500]
     point_colors = rainbow(len(nvalues))
     total_plot = Graphics()
     defint_plot = Graphics()
     for (nidx, nval) in enumerate(nvalues): 
          print "  => N: %d" % nval
          defint = compute_2d_integral(nval, init_vectors)
          beta_ticks = list(np.arange(0.025, 1.025, 0.025))
          #print find_root(defint(alpha) == 0, 4.5, 6)
          defint_plot += plot(defint(alpha), (-5, 5), legend_label = "N = %d" % nval, 
                              legend_color = point_colors[nidx], 
                              rgbcolor = point_colors[nidx])
          alpha_values = map(lambda beta: find_root(defint(alpha) == beta * nval, 0.0001, 4), beta_ticks)
          plot_points = zip(beta_ticks, alpha_values)
          
          print plot_points, "\n"
          xtick_formatter = ["" if n > 0 and n < len(beta_ticks)-1
                            else "%g" % beta_ticks[n] for n in range(0, len(beta_ticks))]
          pplot = points(plot_points, pointsize = 6, 
                         legend_label = "N = %d" % nval, 
                         legend_color = point_colors[nidx], 
                         rgbcolor = point_colors[nidx], 
                         axes_labels = ["$\\beta$", "$\\alpha$"], 
                         gridlines = True, ticks = [beta_ticks, None], 
                         tick_formatter = [xtick_formatter, None], 
                         title = "Looking for Hidden Signals (%s): $\\int_0^1\\int_0^1\\sum_{1 \\leq k \\leq N} \\Re[e^{2\\pi\\imath\\alpha(mx+ny)}] dxdy = \\beta \\cdot N$" % plot_id)
          total_plot += pplot
          total_plot.show()
          defint_plot.show()
     ##
     total_plot.show()
##

def compute_fourier_transform(init_vectors): 
     nvalues = [25, 50, 75, 100, 150, 250, 500, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2500]
     #nvalues = [500, 1000, 1500]
     point_colors = rainbow(len(nvalues))
     total_plot = Graphics()
     slopes, xvals, yvals, scales, ints = [], [], [], [], []
     last_ftsum = lambda omega: 1
     for (nidx, nval) in enumerate(nvalues): 
          ulam_set = compute_ulam_set(nval, init_vectors)
          print "ULAM SET: \n" #, ulam_set, "\n"
          alpha, x, y, omega = var('alpha x y omega')
          ftfunc = lambda m, n: -1 / 4.0 / m / n * (abs(m-omega) + \
                                abs(n - omega) - abs(m + n - omega) \
                                - 2 * abs(omega) + abs(m + omega) + \
                                abs(n + omega) - abs(m + n + omega))
          #ftfunc = lambda m, n: -1 / 4.0 / m / n * (\
          #                      - 1 * abs(omega) + abs(m + omega) + \
          #                      abs(n + omega) - abs(m + n + omega))
          ftfunc2 = lambda m, n: 1 / 4.0 / m / n * (-abs(m-omega) - \
                                abs(n - omega) + abs(m + n - omega) \
                                + abs(m + omega) + \
                                abs(n + omega) - abs(m + n + omega))
          ftfunc3 = lambda m, n: real(-(1.0/(8.0 * m * n)) * (abs(-m + omega) + abs(m + omega) + abs(-n + omega) - 
                                      abs(-m - n + omega) + abs(n + omega) - abs(m + n + omega) - 
                                      i * (-m + omega)**2 * pi * sign(m - omega) - 
                                      i * (-n + omega)**2 * pi * sign(n - omega) + 
                                      i * (-m - n + omega)**2 * pi * sign(m + n - omega) - 
                                      2 * (abs(omega) + i * omega**2 * pi * sign(omega)) + 
                                      i * (m + omega)**2 * pi * sign(m + omega) + 
                                      i * (n + omega)**2 * pi * sign(n + omega) - 
                                      i * (m + n + omega)**2 * pi * sign(m + n + omega)))
          ftfunc4 = lambda m, n: imag(-(1.0/(8.0 * m * n)) * (abs(-m + omega) + abs(m + omega) + abs(-n + omega) - 
                                      abs(-m - n + omega) + abs(n + omega) - abs(m + n + omega) - 
                                      i * (-m + omega)**2 * pi * sign(m - omega) - 
                                      i * (-n + omega)**2 * pi * sign(n - omega) + 
                                      i * (-m - n + omega)**2 * pi * sign(m + n - omega) - 
                                      2 * (abs(omega) + i * omega**2 * pi * sign(omega)) + 
                                      i * (m + omega)**2 * pi * sign(m + omega) + 
                                      i * (n + omega)**2 * pi * sign(n + omega) - 
                                      i * (m + n + omega)**2 * pi * sign(m + n + omega)))
          ftfunc5 = lambda m, n, omega: omega * (m+n)- 2*float(pi)*floor(omega*(m+n) / 2.0 / float(pi))
          ftsum = sum(map(lambda (m, n): ftfunc(m, n) / float(nval * log(nval)), ulam_set))
          #integr = integrate(ftsum, omega, 0, 200) 
          #print integr
          #ints += [(nval, integr)]
          #ftsum2 = ftsum / integr
          ftsum3 = lambda omega: ftsum(omega + 4.745154852900353 * log(nval))
          #print "CRITICAL POINTS: ", solve(derivative(ftsum(omega), omega)==0, omega)
          #print "CRITICAL POINTS II: ", solve(derivative(derivative(ftsum(omega), omega), omega)==0, omega)
          ftplot = plot(ftsum(omega), (-50, 50), legend_label = "N = %d" % nval, 
                        rgbcolor = point_colors[nidx], 
                        title = "S_N-like")
          
          print "omega = 2.4"
          omega = 2.4
          ftset = map(lambda (m, n): ftfunc5(m, n, omega), ulam_set)
          ftplot = histogram(ftset, bins = 250, alpha = 0.5, title = "N = %d, omega = %03g" % (nval, omega), normed=True)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 1.0
          ftset = map(lambda (m, n): cos(omega*(m+n)), ulam_set)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval
          #ftplot.show()
          
          omega = 4.8
          ftset = map(lambda (m, n): ftfunc5(m, n, omega), ulam_set)
          ftplot = histogram(ftset, bins = 250, alpha = 0.5, title = "N = %d, omega = %03g" % (nval, omega), normed=True)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 2.0
          ftset = map(lambda (m, n): cos(omega*(m+n)), ulam_set)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval
          #ftplot.show()
          
          omega = 7.2
          ftset = map(lambda (m, n): ftfunc5(m, n, omega), ulam_set)
          ftplot = histogram(ftset, bins = 250, alpha = 0.5, title = "N = %d, omega = %03g" % (nval, omega), normed=True)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 3.0
          ftset = map(lambda (m, n): cos(omega*(m+n)), ulam_set)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval
          #ftplot.show()
          
          omega = 9.6
          ftset = map(lambda (m, n): ftfunc5(m, n, omega), ulam_set)
          ftplot = histogram(ftset, bins = 250, alpha = 0.5, title = "N = %d, omega = %03g" % (nval, omega), normed=True)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 4.0
          ftset = map(lambda (m, n): cos(omega*(m+n)), ulam_set)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 4.0
          #ftplot.show()
          
          omega = 12.0
          ftset = map(lambda (m, n): ftfunc5(m, n, omega), ulam_set)
          ftplot = histogram(ftset, bins = 250, alpha = 0.5, title = "N = %d, omega = %03g" % (nval, omega), normed=True)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 5.0
          ftset = map(lambda (m, n): cos(omega*(m+n)), ulam_set)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval
          #ftplot.show()
          
          omega = 4.7775
          ftset = map(lambda (m, n): ftfunc5(m, n, omega), ulam_set)
          ftplot = histogram(ftset, bins = 250, alpha = 0.5, title = "N = %d, omega = %03g" % (nval, omega), normed=True)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval / 6.0
          ftset = map(lambda (m, n): cos(omega*(m+n)), ulam_set)
          print "OMEGA(%g), SUM / N: ", sum(ftset), sum(ftset) / nval
          ftplot.show()
          
          #scales += [(nval, sum(ftset) / 6.0 / float(1.094655622762944 * exp(0.1727548405121491 * nval / 100)))]
          scales += [(nval, sum(ftset) / nval)]
          print "SCALES: ", scales, "\n"
          #plotfunc = lambda x: exp((-(x + 4.745154852900353 * log(nval))^2 )+x)
          #ftplot += plot(plotfunc, (-50, 50), legend_label = "N = %d" % nval, 
          #              rgbcolor = point_colors[nidx], 
          #              title = "Fourier Transforms")
          
          #rhsmax = find_local_maximum(ftsum, golden_ratio + 1.5, 75)
          #(x0, y0), (x1, y1) = (float(golden_ratio + 1), ftsum(float(golden_ratio + 1))), \
          #                     (rhsmax[1], ftsum(rhsmax[1]))
          #print "SLOPE: ", (y1-y0) / (x1-x0), "\n"
          #slopes += [(nval, (y1-y0) / (x1-x0))]
          #print "N=% 4d: " % nval, find_local_maximum(ftsum,-2.7,-1.6), find_local_minimum(ftsum,-2.7,-1.6)
          #print "N=% 4d: " % nval, find_local_maximum(ftsum,1.6,2.7), find_local_minimum(ftsum,1.6,2.7)
          #print "N=% 4d: " % nval, vector(find_local_maximum(ftsum,1.6,2.7))/nval, vector(find_local_minimum(ftsum,1.6,2.7))/nval
          #print "N=% 4d: " % nval, vector(find_local_maximum(ftsum,-75,-6))/nval, vector(find_local_minimum(ftsum,-75,-6)) / nval
          #print "N=% 4d: " % nval, vector(find_local_maximum(ftsum,6,75)) / nval, vector(find_local_minimum(ftsum,6,75)) / nval
          #print ""
          total_plot += ftplot
          #total_plot.show()
          #print float(last_ftsum(0)), float(ftsum3(0))
          #scales += [(nval, float(ftsum3(0)) / float(last_ftsum(0)))]
          #scales += [(nval, float(ftsum3(0)))]
          #last_ftsum = ftsum3(omega)
          #print "SCALES: ", scales, "\n"
          #print "INTS: ", ints, "\n"
     ## 
     #print "INTS: ", ints, "\n"
     #total_plot += list_plot(ints, rgbcolor='red')
     total_plot = list_plot(scales, rgbcolor='blue', plotjoined=True)
     total_plot.show()
##

initial_vector_configs = [ ## examples from Jayadev's talk and in the article: 
          [V(1, golden_ratio), V(0, 1)], 
          [V(1, golden_ratio), V(golden_ratio, 1)], 
          [V(1, golden_ratio), V(1, 0)], 
          #[V(1, 0), V(0, 1)], 
          #[V(9, 0), V(0, 9), V(1, 13)], 
          #[V(2, 5), V(3, 1)], 
          #[V(1, 0), V(2, 0), V(0, 1)], 
          #[V(2, 0), V(0, 1), V(3, 1)], 
          #[V(1, 0), V(0, 1), V(2, 3)], 
          #[V(3, 0), V(0, 1), V(1, 1)], 
          #[V(1, 0), V(2, 0), V(0, 1)], 
          #[V(2, 0), V(3, 0), V(0, 1)], 
          #[V(1, 0), V(0, 1), V(6, 4)], 
          #[V(1, 0), V(0, 1), V(10, 9)], 
          #[V(1, 0), V(0, 1), V(10, 3)], 
          #[V(1, 3), V(3, 4)], 
          #[V(1, 0), V(1, 1)]
     ] 
Nvalue = 20
init_vectors = initial_vector_configs[0]


     
