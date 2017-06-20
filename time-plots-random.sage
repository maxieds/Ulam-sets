from sys import argv
from oset import oset ## Source distribution: https://pypi.python.org/pypi/oset/
from random import uniform, seed
from datetime import datetime
from sage.all import *
from sage.plot.colors import rainbow
from ulam_sets_modified import compute_ulam_set
from Means import MeanStats

X = lambda v: v[0]
Y = lambda v: v[1]

def edist(v, w): 
     vdim = 1 if not isinstance(v, tuple) else len(v)
     [vsum, wsum] = map(sum, [v, w])
     return vsum**2 + wsum**2
##

def generate_random_vectors(vector_1d = False, intmin = 0, intmax = 1, numvecs = 2, exclude_zero = True): 
     random.seed(datetime.now())
     random_func = lambda: random.uniform(intmin, intmax)
     draw_vector_func = lambda: random_func() if vector_1d else (random_func(), random_func())
     is_zero = lambda vec: vec == 0  if not isinstance(v, tuple) else vec == (0, 0)
     init_vectors, rvec_count, prevvec = [], 0, None
     while rvec_count < numvecs: 
          rvec = tuple(map(float, draw_vector_func()))
          if not is_zero(rvec) or not exclude_zero and (rvec_count == 0 or rvec != prevvec): 
               init_vectors += [rvec]
               rvec_count += 1
               prevvec = rvec
          ##
     ##
     return init_vectors
## 

def generate_continuous_histogram_plot(datapts, bins = 250, cdf = False, normed = True, rgb = 'blue'): 
     histpts, bin_edges = np.histogram(datapts, bins = bins, density = normed, 
                                       cumulative = cdf, label = None)
     bincenters = np.diff(bin_edges)
     plotpts = zip(list(bincenters), list(histpts))
     rplot = list_plot(plotpts, plotjoined = True, rgbcolor = rgb, legend_label = label, 
                       fill = axis, fillcolor = rgb, fillalpha = 0.5, thickness = 2)
     rplot.set_legend_options(title = "Statistics for Times", fancybox = True, 
                              font_size = 'large', ncol = 2, back_color = '#c2c2c2')
     return rplot
##

def compute_nthline_points(n, ulam_set, v0, v1, v2): ## TODO: Check this function for correctness 
     print "v0, v1, v2: ", v0, v1, v2
     slope_func = lambda (x0, y0), (x1, y1): float((y1-y0) / (x1-x0))
     upper_line = lambda x: x if v0[0]-v2[0] == 0 else slope_func(v0, v2) * (x - X(v0)) + Y(v0)
     lower_line = lambda x: v1[1] if v1[0]-v2[0] == 0 or v1[1]-v2[1] == 0 \
                                  else slope_func(v1, v2) * (x - X(v1)) + Y(v1)
     upperhyp, lowerhyp = n * edist(v0, v2), n * edist(v1, v2)
     upper_linex = lambda x: 0 if v0[0]-v2[0] == 0 else float(sqrt(upperhyp**2 - upper_line(x)**2) + X(v2))
     lower_linex = lambda x: x if v1[1]-v2[1] == 0 else float(sqrt(lowerhyp**2 - lower_line(x)**2) + X(v2))
     
     #print eval_on_operands(upper_line)(z)
     #var('z')
     #print lower_line(z), upper_line(z)
     
     upper_line_vec = lambda n: (upper_linex(n), upper_line(n))
     lower_line_vec = lambda n: (lower_linex(n), lower_line(n))
     sv0, sv1 = upper_line_vec(n), lower_line_vec(n)
     
     if n == 0: 
          return [v2]
     elif n % 2 == 1: # n odd: only the boundary vectors here
          return [sv0, sv1]
     else: 
          numpts = (n + 2.0) / 2.0
          segment_line = lambda x: slope_func(sv0, sv1) * (x - X(sv0)) + Y(sv0)
          return [segment_line(x) for x in range(0, numpts)]
     ##
##

def compute_nthline_full(n, ulam_set, v0, v1, v2): 
     searchpts = compute_nthline_points(n, ulam_set, v0, v1, v2)
     print "n, searchpts: ", n, searchpts
     nthline_full = True
     for sp in searchpts: 
          print "filter: ", sp, filter(lambda (ntimes, vec): vec == sp, ulam_set)
          if len(filter(lambda (ntimes, vec): vec == sp, ulam_set)) != 1: 
               nthline_full = False
               break
          ##
     ##
     print "nthline_full: ", nthline_full
     return nthline_full
##

def compute_nthline_times(ulam_set, v0, v1, v2): 
     nsteps, timeslst = 1, []
     while True: 
          if not compute_nthline_full(nsteps, ulam_set, v0, v1, v2): 
               break
          ##
          searchpts = compute_nthline_points(nsteps, ulam_set, v0, v1, v2)
          nthline_indices = []
          for sp in searchpts: 
               (ntime, vec) = filter(lambda (nidx, vec): vec == sp, ulam_set)[0]
               nthline_indices += [ntime]
          ##
          timeslst += nthline_indices
     ##
     return timeslst
## 

## Usage: <prog-name> METHOD coord1 coord2 [coord3 coord4] NSTEPS
## Example: <prog-name> userdef-1d 1 1 250
## Example: <prog-name> userdef-2d 1 0 0 1 250
## Example: <prog-name> random-1d 250
## Example: <prog-name> random-2d 250
if __name__ == '__main__':
     
     init_vector_method = argv[1]
     init_vector_coords = map(float, argv[2:-1]) + [int(argv[-1])]
     init_vectors, nsteps = [], 0
     if init_vector_method == "userdef-1d": 
          init_vectors += [init_vector_coords[0], init_vector_coords[1]]
          nsteps = init_vector_coords[2]
     elif init_vector_method == "userdef-2d": 
          v0 = (init_vector_coords[0], init_vector_coords[1])
          v1 = (init_vector_coords[2], init_vector_coords[3])
          init_vectors += [v0, v1]
          nsteps = init_vector_coords[4]
     elif init_vector_method == "random-1d": 
          init_vectors = generate_random_vectors(vector_1d = True, numvecs = 2)
          nsteps = argv[2]
     elif init_vector_method == "random-2d" or init_vector_method == "none": 
          init_vectors = generate_random_vectors(vector_1d = False, numvecs = 2)
          nsteps = argv[2]
     else: 
          raise ValueError("Unsupported initial vector method.")
     ##
     
     v0, v1 = init_vectors
     wedge_angle = abs(math.atan2(Y(v0), X(v0)) - math.atan2(Y(v1), X(v1)))
     ulam_set = compute_ulam_set(nsteps, init_vectors = [v0, v1], return_indices = True)
     #sys.exit(0)
     nthline_times = compute_nthline_times(ulam_set, v0, v1, v2 = ulam_set[2][1])
     nthline_stats = map(lambda nthlst: MeanStats.compute_all_statistics(nthlst), nthline_times)
     print "nthline_times: ", nthline_times
     print "nthline_stats: ", nthline_stats
     statdesc = map(lambda (sdesc, stat): sdesc, nthline_stats[0])
     point_colors = rainbow(len(statdesc))
     total_plot = Graphics()
     for (sidx, sdesc) in enumerate(statdesc): 
          getnthline_stats = nthline_stats[:][sidx]
          plotn = generate_continuous_histogram_plot(getnthline_stats, label = sdesc, rgb = point_colors[sidx])
          total_plot += plotn
     ##
     
     plot_title_vlst = map(lambda (x, y): r'$[^{% 3g}_{% 3g}]$' % (x, y), [v1, v2])
     total_plot_title = "ULAM SET $\\equiv$ \{" + ', '.join(plot_title_vlst) + ", \ldots\}"
     image_path = "random-ulam-set-" + str(datetime.now())
     total_plot.save(image_path, title = total_plot_title, frame = True, 
                     typeset = 'latex', axes_labels = ("N", "pdf/cdf(N)"))
     total_plot.show(title = total_plot_title, frame = True, 
                     typeset = 'latex', axes_labels = ("N", "pdf/cdf(N)"))
##
     
     
     
     
