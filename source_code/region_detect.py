
import numpy as np
import icesatUtils as it


def calc_overlap(x1, x2, dx, dx_overlap=0.0, dx2=None):
    """
    Find the inner overlaps between arrays x1 and x2
    Input:
        x1 - 1d array
        x2 - 1d array
        dx - change in x assigned to each point
        dx_overlap - addition on either end of each region
        dx2 - optional dx for x2, if sampling rates are different
                i.e. overlaps in 03 and 08

    Output
        return regions_overlap
        regions of start/end ttg

    Example:
        import region_detect as rd
        # ... load in 03 data
        dt0 = 1.0 / 10000 # change in time for ATL03
        regions_overlap = rd.calc_overlap(t1, t2, 10*dt0, dt0)

    """

    dt1 = dx
    regions2, _ = find_region(x2, dx, dx_overlap)
    regions2 = combine_region(regions2)

    if type(dx2) == type(None):
        dx2 = dx
    regions1, _ = find_region(x1, dx2, dx_overlap)
    regions1 = combine_region(regions1)

    regions_overlap = []
    for reg2 in regions2:
        for reg1 in regions1:
            reg = inner_region(reg2, reg1)
            if reg != []:
                regions_overlap.append(reg)

    return regions_overlap

def filt_reg_x(x, regions, equal=False):
    """
    Given a list of regions, this
    filters the array x for only values
    in those regions.

    Input:
        x - numpy array of values
        regions - list of 1d start/end values
                i.e., [[t0,t1], [t2,t3], ...]
        equal - boolean to include end-points in 1d
                regions or not

    Output:
        the filtered x values and corresponding
        boolean indexing, which may concievably
        be applied to some corresponding y values
        return x_filt, jr

    Example:
        import numpy as np
        import region_detect as rd
        x = np.linspace(0,10)
        y = np.linspace(10,20)
        regions = [[0,3], [6,7]]
        x_filt, jr = rd.filt_reg_x(x, regions)

        y_filt = y[jr]

        # x_filt contains values only within (0,3)
        # and (6,7)

    """

    if type(x) != np.ndarray:
        x = np.array(x)

    if type(regions) != list:
        regions = np.ndarray.tolist(regions)

    jr = np.zeros(len(x)).astype(bool)
    for reg in regions:
        if equal:
            b = (x >= reg[0]) & (x <= reg[1])
        else:
            b = (x > reg[0]) & (x < reg[1])
        jr[b] = True
    x_filt = x[jr]
    return x_filt, jr

def filt_reg(x, y, regions, equal=False):

    """
    Wrapper for filt_reg_x(), additionally
    returning y_filt

    See filt_reg_x for input/output/example details.
    """

    if type(x) != np.ndarray:
        x = np.array(x)
    if type(y) != np.ndarray:
        y = np.array(y)
    if type(regions) != list:
        regions = np.ndarray.tolist(regions)

    x_filt, jr = filt_reg_x(x, regions, equal)
    y_filt = y[jr]

    return x_filt, y_filt, jr


def density_detect(x, z, r_min=1e-3, r_max=25.0, loading_bar=False, num_cmp=None):

    """
    This function creates a density map out of point data,
    using nearest-neighbors-like operation. This function is meant
    to be used as a quick-look for a density profile, or for 
    intially obtaining training data for machine-learning 
    applications.

    The method is similar to a clustering algorithm,
    which goes to each point, looks at a small region
    of 2D space, and assigns a "density" to each point
    based on the inverse of the distance of other points.

    Input:
        x - array of domain
        z - array of range
        r_min - minimum radius to calc density
        r_max - maximum radius to calc density
        loading_bar - boolean to activate tqdm loading
        num_cmp - option (int) to cap how many points
                    can be used in density calc
    
    Output:
        Returns density, rho, and number of points used
        in each density calc

        return rho, num_pts

    Example:
        import numpy as np
        import icesatReader as ir
        import icesatPlot as ip
        import region_detect as rd
        file_03 = '/bigtex_data/data/release/003/ATL03_r003' + \
                    '/Alaska/ATL03_20200118153211_03500605_003_01.h5'
        file_08 = '/bigtex_data/data/release/003/ATL08_r003' + \
                    '/Alaska/ATL08_20200118153211_03500605_003_01.h5'
        gt = 'gt1r'
        atl03 = ir.get_atl03_struct(file_03, gt, file_08)

        df_03 = atl03.df

        # remove duplicate point/times
        t = np.array(df_03.delta_time)
        x = np.array(df_03.alongtrack)
        z = np.array(df_03.h_ph)
        t_uq, index = np.unique(t, return_index=True)
        x, z = x[index], z[index]

        # filter for a small region
        x, z, _ = rd.filt_reg(x, z, [[150000, 200000]])

        rho, num_pts = rd.density_detect(x, z, r_max=5.0, loading_bar=True)

        fig, ax = ip.make_fig()
        ax.scatter(x, z, s=10, c=rho)
        fig.show()

    Caveats:
        Does not handle effects on ends of data;
        that is, density will appear less on the
        edges since data does not exist further.

        x and z must have the same units for
        best results.

    """
    limit_cmp = False
    if type(num_cmp) != type(None):
        limit_cmp = True
        num_cmp = int(num_cmp)

    n = len(x)
    rho = np.zeros(n)
    num_pts = np.zeros(n)

    from tqdm import tqdm
    iterator = range(n)
    if loading_bar:
        from tqdm import tqdm
        iterator = tqdm(range(n))

    for j in iterator:

        x1, x2, x3 = x[j]-r_max, x[j], x[j]+r_max

        dj = 0
        while True:
            if j-dj > 0:
                if not (x[j-dj] > x1):
                    dj += 1
                    break
            else:
                break
            dj += 1
        dj_lower = dj

        dj = 0
        while True:
            if j+dj < n:
                if not (x[j+dj] < x3):
                    dj -= 1
                    break
            else:
                break
            dj += 1
        dj_upper = dj


        x_reg, z_reg = x[j-dj_lower:j+dj_upper+1], z[j-dj_lower:j+dj_upper+1]
        # if abs(x_reg[0] - x1) > dr:
        #   continue
        # elif abs(x_reg[1] - x3) > dr:
        #   continue

        if len(x_reg) == 0:
            continue

        dx_reg = abs(x_reg - x[j])
        dz_reg = abs(z_reg - z[j])

        r_reg = np.sqrt(dx_reg**2 + dz_reg**2)
        r_reg = r_reg[(r_reg > r_min) & (r_reg < r_max)]

        if len(r_reg) == 0:
            continue
        
        # up to num_cmp points considered
        if limit_cmp:
            r_reg = np.sort(r_reg)
            nr = min(len(r_reg), num_cmp)
            r_reg = r_reg[:nr]

        num_pts[j] = len(r_reg)
        k_reg = 1.0 / r_reg

        # rho[j] = np.mean(k_reg)
        rho[j] = sum(k_reg)
        # rho[j] = np.log10(sum(k_reg))


    return rho, num_pts


def path_detect(t, y, w, dt_max, dy_max, dt_min=-np.inf, dy_min=-np.inf, 
                w_id=None, flip_filter=False, rtn_id=False, loading_bar=False, debug=0):

    """
    Conceptually, this function draws a path using the closest 
    new point to end of the current path. If points are not
    found as part of this path, a new path is made.

    The default return is all paths as one object, but
    rtn_id = True will enable the return of all different paths
    detected.

    This function is complete, but limited relative to the 
    filtered_centroids.py case. This function returns j_valid,
    which is a list of indices for all paths found. Effectively,
    this would find elements of t that are -not- noise, instead
    of finding particular paths. However, this function
    could be adapted to finding particular paths (per object),
    such as a path for star 1 and another for star 2 and so on.

    Input:
    t - time-series, that has no negative entries
    y - some data corresponding to time-tags
    w
        - local window size
    dt_max
        - limit to dt between path points
            (i.e. noise may be higher)
    dy_max
        - limit to dy between path points
            (i.e. noise may be higher)
    dt_min/dy_min
        - extra limit on how to choose next point
    w_id
        - window size required to identify these points as
            a path
        - default is w/2
    flip_filter
        - look for what is farthest, rather than closest
    rtn_id
        - return the id of each path

    Output:
    j_valid
        - indices where paths exist, or where signal of
            some kind exists

    Caveats:
        The local windowing mechanism keeps the algorithm roughly
        linear, but j_valid is continuously increasing, at a rapid
        pace, especially when w is large. Advise using either 
        sections of data, or less than 1e6 points. I think:
            if w = 100 and len(t) = 1e6, then it will take T time
            if w = 200 and len(t) = 1e6, then it will take 2*T time
            if w = 100 and len(t) = 1e5, then idk what time it will take

    Future Updates:
        Allow for user modification of the detection filter, as some
        arbitrary function for dt and dy. That way, you could make 
        it work for the magnitude of dt and dy (like filter_centroids.py)
        or for negative dt.

    Example:
        # prep ground points to be separated by SVM
        t_seg_ground = t_seg[c_seg == 1]
        h_seg_ground = h_seg[c_seg == 1]
        dt = 1.0 / 10000.0
        j_path, j_id = rd.path_detect(t_seg_ground, h_seg_ground, w=50, dt_max=5*dt, dy_max=0.3, \
                                 dt_min=-np.inf, dy_min=-np.inf, \
                                 w_id=None, flip_filter=False, rtn_id=True, loading_bar=False, debug=3)

    """

    class Object():
        def __init__(self, tau0, y0, j0, object_id):
            self.t = [tau0]
            self.y = [y0]
            self.j = [j0]
            self.object_id = object_id
            self.init = 0
            self.tau = 0.0

        def limit(self):
            if len(self.t) > w:
                del self.t[0]
                del self.y[0]
                del self.j[0]

        def remove_first(self):
            if len(self.t) > 0:
                del self.t[0]
                del self.y[0]
                del self.j[0]


    obj_list = []
    if w_id == None:
        w_id = int(w/2)

    object_id = 0
    
    n = len(t)
    j_path = np.zeros(n).astype(bool)
    j_id = {}

    iterator = range(n)
    if loading_bar:
        from tqdm import tqdm
        iterator = tqdm(range(n))

    # for j in range(n):
    for j in iterator:

        # m = -1
        # if len(obj_list) > 0:
        #   m = min([len(obj_list[i].t) for i in range(len(obj_list))])
        # out_str = '%d, %d' % (len(obj_list), m)
        # sys.stdout.write(out_str)
        # sys.stdout.flush()
        # backspace = ''
        # for k in range(len(out_str)):
        #   backspace += '\b'
        # sys.stdout.write(backspace)

        tau0, y0 = t[j], y[j]
        if len(obj_list) == 0:
            obj_list.append(Object(tau0, y0, j, object_id))
        else:
            t_new, y_new = t[j], y[j]
            possible_match = []
            for i, obj in enumerate(obj_list):
                # num_pts = len(obj.t)
                # dj = min(num_pts, 3)
                # t_last_arr, y_last_arr = obj.t[-dj:], obj.y[-dj:]

                # for j_rel in range(dj):
                #   t_last, y_last = t_last_arr[j_rel], y_last_arr[j_rel]
                #   dt = t_new - t_last
                #   if dt > dt_max:
                #       b = 1
                #   dy_abs = abs(y_last - y_new)
                #   if dy_abs > dy_max:
                #       b = 1

                #   if b:
                #       continue

                #   possible_match.append([dy_abs, i])


                t_last, y_last = obj.t[-1], obj.y[-1]

                b = 0
                # detection filter
                ####################
                if not flip_filter:
                    dt = t_new - t_last
                    if dt > dt_max:
                    # if not (dt_min < dt < dt_max):
                        b = 1
                    dy_abs = abs(y_last - y_new)
                    if dy_abs > dy_max:
                    # if not (dy_min < dy_abs < dy_max):
                        b = 1
                else:
                    dt = t_new - t_last
                    if dt < dt_max:
                    # if (dt_min < dt < dt_max):
                        b = 1
                    dy_abs = abs(y_last - y_new)
                    if dy_abs < dy_max:
                    # if (dy_min < dy_abs < dy_max):
                        b = 1
                ####################
                if b:
                    continue

                # match by abs(dh) elevation change within
                # valid time change dt_max
                possible_match.append([dy_abs, i])


            if len(possible_match) > 0:
                # there exists at least one possible match
                possible_match = sorted(possible_match)
                i0 = possible_match[0][1]
                obj = obj_list[i0]

                obj.t.append(t_new)
                obj.y.append(y_new)
                obj.j.append(j)

                # remove first element if the object has
                # reached its window length
                # obj.limit()

                num_pts = len(obj.t)
                if num_pts > w_id:
                    # if obj.init == 0:
                    #   obj.tau = obj.t[0]
                    #   obj.init = 1

                    # else:
                    #   j0 = np.argmin(abs(np.array(obj.t) - obj.tau))
                    #   for jj in range(j0+1,num_pts):
                    #       j_path.append(jj)

                    # j_id0 = j_id[obj.object_id]

                    if rtn_id:
                        if not (obj.object_id in j_id):
                            j_id[obj.object_id] = []

                    for jj in obj.j:
                        # j_path.append(jj)
                        j_path[jj] = True
                        if rtn_id:
                            if len(j_id[obj.object_id]) > 0:
                                j_last = j_id[obj.object_id][-1]
                                if jj > j_last:
                                    j_id[obj.object_id].append(jj)
                            else:
                                j_id[obj.object_id] = [jj]


            else:
                # create new object, point is too far from any objects
                object_id += 1
                obj_list.append(Object(t_new, y_new, j, object_id))

        # remove the first element off each object,
        # whether they match or not
        #   noise will be removed
        i = 0
        while i < len(obj_list):
            obj = obj_list[i]
            if t[j] - obj.t[-1] > 2*dt_max: # 2*dt, safe side
                obj.remove_first()
                if len(obj.t) == 0:
                    del obj_list[i]
                    i -= 1
            i += 1

        # for i, obj in enumerate(obj_list):
        #   if len(obj.t) > w_id:
        #       for jj in obj.j:
        #           j_path.append(jj)

    # if debug:
    #   print('len(j_path)', len(j_path))
    # j_path = np.unique(j_path)

    if rtn_id:
        return j_path, j_id
    else:
        return j_path



def est_dx(x, rtn_s=False, IQR_factor=1.5, debug=0):

    """
    Estimates the change in x, while removing outliers in that change.
    rtn_s
        - returns both dx, and standard deviation of dx
        - return dx, dx_s

    IQR_factor
        - default is 1.5, for identifying outliers
    
    Output:
        if rtn_s:
            return dx, dx_s
        else:
            return dx

    """

    dx_range = np.diff(x)
    if len(dx_range) == 0:
        if debug:
            print('warning: len(dx_range) == 0')
        return 0


    Q1,Q2,Q3,IQR = it.get_outlier_data(dx_range)
    dx_range_red = []

    # assumes dx_range to have all positive values
    # abs(dx_range[jj]) is used in filter loop, so it's
    # consistent if dx_range has negative values
    m1, m2 = 0, abs(Q1-IQR_factor*IQR) + abs(Q3+IQR_factor*IQR)
    # if sign == '+':
    #   # if x is known to always be increasing
    #   m1, m2 = 0, abs(Q1-IQR_factor*IQR) + abs(Q3+IQR_factor*IQR)
    # else:
    #   m1, m2 = Q1-IQR_factor*IQR, Q3+IQR_factor*IQR
    debug_index = []
    for jj in range(len(dx_range)):
        if m1 < abs(dx_range[jj]) < m2:
            dx_range_red.append(dx_range[jj])
            debug_index.append(jj)

    if debug > 1:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(range(len(dx_range)), dx_range)
        ax.plot([debug_index[0], debug_index[-1]], [m1, m1], 'r--')
        ax.plot([debug_index[0], debug_index[-1]], [m2, m2], 'r--')

    if len(dx_range_red) == 0:
        if debug:
            print('warning: len(dx_range_red) == 0')
        dx = np.nanmean(dx_range)
        dx_s = np.nanstd(dx_range)
    else:
        dx = np.nanmean(dx_range_red)
        dx_s = np.nanstd(dx_range_red)

    if debug > 1:
        ax.plot([debug_index[0], debug_index[-1]], [dx, dx], '--', color='C1')
        ax.plot([debug_index[0], debug_index[-1]], [dx + 3*dx_s, dx + 3*dx_s], '--', color='C1')
        fig.show()

    if rtn_s:
        return dx, dx_s
    else:
        return dx


def overlap(reg1, reg2, equal=True):

    """
    reg1 and reg2 are lists, such as
        reg1 = [t1, t2]
        reg2 = [t3, t4]

    equal=True means that if
      reg1 = [10,20]
      reg2 = [20,40],
    output will be 1 since 20's are the same
    """
    if type(reg1) != list:
        reg1 = list(reg1)
    if type(reg2) != list:
        reg2 = list(reg2)

    reg = sorted([reg1, reg2])
    reg1_t, reg2_t = reg[0], reg[1]

    b = 0
    if equal:
        if reg1_t[1] < reg2_t[0]:
            return 0
        else:
            return 1
    else:
        if reg1_t[1] <= reg2_t[0]:
            return 0
        else:
            return 1

def inner_region(reg1, reg2, equal=True):
    if type(reg1) != list:
        reg1 = list(reg1)
    if type(reg2) != list:
        reg2 = list(reg2)

    if overlap(reg1, reg2, equal):
        total_reg = sorted([reg1[0], reg1[1], reg2[0], reg2[1]])
        return [total_reg[1], total_reg[2]]
    else:
        return []

def eq_check(reg):
    if reg[1] - reg[0] != 0:
        return reg
    else:
        return []


def mask_region(A, B, priority, equal=True):
    
    """
    This takes in list regions A and B, such as
        A = [t1, t2]
        B = [t3, t4],
    and masks one region over the other, depending
    on the priority. If priority == 1, then A stays
    the same, and B gets cut. If priority == 2, then
    B stays the same and A gets cut to fit B.

    Input:
        A - first region, such as [t1, t2]
        B - second region, such as [t3, t4]
        priority - int [1 or 2] which defines which
                    region stays and the other gets
                    masked

    Output:
        basically, the new masked regions
        return A_new, B_new

    """

    if type(A) != list:
        A = list(A)
    if type(B) != list:
        B = list(B)

    if not priority in [1,2]:
        raise ValueError('priority must be either 1 or 2')

    if A == B:
        if priority == 1:
            return A, []
        elif priority == 2:
            return [], A

    if overlap(A, B, equal):
        # overlap exists
        b = 0
        if equal:
            if B[0] <= A[0] and B[1] >= A[1]:
                b = 1
            elif A[0] <= B[0] and A[1] >= B[1]:
                b = 2
            else:
                b = 3
        else:
            if B[0] < A[0] and B[1] > A[1]:
                b = 1
            elif A[0] < B[0] and A[1] > B[1]:
                b = 2
            else:
                b = 3

        # if B[0] <= A[0] and B[1] >= A[1]:
        if b == 1:
            # full overlap, superset B
            if priority == 1:
                A_new_eq = eq_check([A[0],A[1]])
                B_new = [[B[0],A[0]], [A[1],B[1]]]

                B_new_eq = B_new
                if eq_check(B_new[0]) == []:
                    B_new_eq = B_new[1]
                    if eq_check(B_new[1]) == []:
                        B_new_eq = []

                if eq_check(B_new[1]) == []:
                    B_new_eq = B_new[0]
                    if eq_check(B_new[0]) == []:
                        B_new_eq = []

                return A_new_eq, B_new_eq

            elif priority == 2:
                A_new_eq = []
                B_new_eq = eq_check([B[0],B[1]])
                return A_new_eq, B_new_eq


        # elif A[0] <= B[0] and A[1] >= B[1]:
        elif b == 2:
            # full overlap, superset A
            if priority == 1:
                A_new_eq = eq_check([A[0],A[1]])
                B_new_eq = []
                return A_new_eq, B_new_eq

            elif priority == 2:
                A_new = [[A[0],B[0]], [B[1],A[1]]]
                B_new_eq = eq_check([B[0],B[1]])

                A_new_eq = A_new
                if eq_check(A_new[0]) == []:
                    A_new_eq = A_new[1]
                    if eq_check(A_new[1]) == []:
                        A_new_eq = []

                if eq_check(A_new[1]) == []:
                    A_new_eq = A_new[0]
                    if eq_check(A_new[0]) == []:
                        A_new_eq = []

                return A_new_eq, B_new_eq

        else:
            # partial overlap
            if A[0] < B[0]:
                if priority == 1:
                    A_new = [A[0],A[1]]
                    B_new = [A[1],B[1]]
                    return eq_check(A_new), eq_check(B_new)
                elif priority == 2:
                    A_new = [A[0],B[0]]
                    B_new = [B[0],B[1]]
                    return eq_check(A_new), eq_check(B_new)

            elif B[0] < A[0]:
                if priority == 1:
                    A_new = [A[0],A[1]]
                    B_new = [B[0],A[0]]
                    return eq_check(A_new), eq_check(B_new)
                elif priority == 2:
                    A_new = [B[1],A[1]]
                    B_new = [B[0],B[1]]
                    return eq_check(A_new), eq_check(B_new)

    else:
        # no overlap
        return [], []


def mask_region_multi(region_ids, regions, P, debug=False):


    """
    This function takes in multiple regions defined by start/end
    times and multiple priorities for each region (can have same
    priority for multiple regions), and outputs the highest
    priority overlap for all regions, i.e.
        if given regions A, B
        A: [10, 20], priority == 1 (P[0] = 1)
        B: [15, 17], priority == 2 (P[1] = 2),
        
        output
        A: [[10, 15], [17, 20]]
        B: [[15, 17]]
        since B has a higher priority than A

    Input:
    1. region_ids:
        - some kind of unique identifier for each region
        - can be number, letter, word, etc, outputs as dictionary keys
    2. regions:
        - multiple regions defined like
            regions = [[10, 20], [15, 17], [30, 40], ...]
            note that second number > first number
    3. P:
        - integer priorities for each region (can be negative)
            P = [1, 1, 3, 2, 3, 1, 5, ...]


    Output:
    1. namelist:
        - dictionary of region_ids with entries as the final regions
        - note that a final region can consist of multiple regions for
            the same region id
        - note further than some regions can be removed if they are
            entirely overlapped by a higher-priority region

    2. regions_all: (mostly for debug)
        - all regions in namelist, sorted such that any overlaps
            would be very visible (there should be no overlaps
            at the end)

    """


    num_regions = len(region_ids)
    err = 0
    if num_regions != len(regions):
        print('error: len(region_ids) != len(regions)')
        err = 1
    if num_regions != len(P):
        print('error: len(region_ids) != len(P)')
        err = 1

    # check that regions is correctly formatted,
    # and sort each region by least-to-greatest
    # just in case
    for i in range(num_regions):
        if len(regions[i]) != 2:
            err = 1
            print('error: len(regions[%d]) != 2' % i)
            # break
        regions[i] = sorted(regions[i])

        if type(P[i]) != int:
            # break
            p = int(P[i])
            if p != P[i]:
                err = 1
                print('error: type(P[%d]) != int' % i)
            P[i] = p

    if err:
        return {}, []

    ############################################
    # # Testing
    # seed = int(sys.argv[1])
    # np.random.seed(seed)
    # num_regions = 75
    # # max_time = 150
    # d_time = 25
    # max_priority = 15

    # P = []
    # regions = []
    # for i in range(num_regions):
    #   while 1:
    #       # a = np.random.randint(max_time)
    #       # b = np.random.randint(max_time)
    #       # t0, tf = i*5, i*5 + d_time
    #       t0, tf = i, i + d_time
    #       # a = np.random.randint(t0,tf)
    #       # b = np.random.randint(t0,tf)
    #       a = np.random.rand(1)[0]*(tf-t0) + t0
    #       b = np.random.rand(1)[0]*(tf-t0) + t0
    #       if abs(a) > abs(b+1):
    #           break

    #   regions.append(sorted([a,b]))
    #   P.append(np.random.randint(max_priority))

    # a0 = 65
    # # region_ids = [str(unichr(a0+i)) for i in range(num_regions)]
    # region_ids = [str(i) for i in range(num_regions)]

    # namelist, regions_all = rd.mask_region_multi(region_ids, regions, P, debug=False)
    ############################################



    # sort by starting time
    regions_p = [[regions[i],P[i],region_ids[i]] for i in range(num_regions)]
    regions_p.sort(key=lambda x: x[0][0])
    regions = [regions_p[i][0] for i in range(num_regions)]
    regions_check = combine_region(regions) # for debug
    P = [regions_p[i][1] for i in range(num_regions)]
    region_ids = [regions_p[i][2] for i in range(num_regions)]


    # increase duplicate entries in P so there is no
    # equal priority
    #   sort P, then increase remaining list if duplicates are found
    #   i.e. given P       = [2, 4, 1, 2, 1, 2, 5, 5, 2]  (after regions time-sorting)
    #              P_s     = [1, 1, 2, 2, 2, 2, 4, 5, 5]  (sort by priority)
    #              P_s_new = [1, 2, 3, 4, 5, 6, 8, 9, 10] (apply algorithm)
    #              P_new   = [3, 8, 1, 4, 2, 5, 9, 10, 6] (reverse sort back to original P)

    P0 = list(np.copy(P))
    # print('P', P)
    if len(np.unique(P)) != len(P):
        # some duplicate entries exist

        P_old = list(np.copy(P))
        P_s_index = np.argsort(P_old)   # save indicies for reversing sort later
        P_s0 = list(np.array(P_old)[P_s_index])
        P_s = list(np.copy(P_s0))
        # print('P_s', P_s)

        P_s_new = []
        dp = 0
        for i in range(num_regions):
            if not (P_s[i] in P_s_new):
                P_s_new.append(P_s[i])
            else:
                P_s = P_s[:i] + list(np.array(P_s[i:]) + 1)
                P_s_new.append(P_s[i])

        # print('P_s_new', P_s_new)

        P_new = np.zeros(num_regions).astype(int)
        for i in range(num_regions):
            P_new[P_s_index[i]] = P_s_new[i]
        P_new = list(P_new)

        # print('P_new', P_new)
        P = P_new

    # sys.exit()

    namelist = {}
    # for i, f in enumerate(files):
    #   namelist[f] = [regions[i]]

    if debug:
        print(regions)
        print(P0)
        print(P)
        print('')


    # loop over all regions and determine how region i is limited
    # by all regions j, where j != i
    #   this makes reg_limit_p, regions that region i cannot be in

    for i, reg1 in enumerate(regions):

        reg_limit_p = []
        for j, reg2 in enumerate(regions):
            if i == j:
                continue

            if overlap(reg1, reg2):
                if P[j] > P[i]:
                    # reg1_new, reg2_new = overlap_priority(reg1, reg2, 2)
                    reg1_new, reg2_new = mask_region(reg1, reg2, 2)
                    if len(reg2_new) > 0:
                        if type(reg2_new[0]) == list:
                            # reg2_new == [[num, num], [num, num]]
                            reg_limit_p.append([reg2_new[0], P[j]])
                            reg_limit_p.append([reg2_new[1], P[j]])
                        else:
                            # reg2_new == [num, num]
                            reg_limit_p.append([reg2_new, P[j]])



        # limit by highest priority first
        reg_limit_p.sort(key=lambda x: x[1], reverse=True)
        reg_limit = [arr[0] for arr in reg_limit_p] # regions that region i cannot be in
        if debug:
            print(i, P[i], reg_limit_p)

        reg1_total = [reg1]
        if len(reg_limit_p) > 0:
            # if limited at all
            k = 0
            while k < len(reg1_total):
                reg1_t = reg1_total[k]
                # if reg1_t == []:
                #   del reg1_total[k]
                #   # k -= 1
                #   # k += 1
                #   continue

                for reg_lim in reg_limit:
                    if overlap(reg1_t, reg_lim):
                        # reg1_t, _ = overlap_priority(reg1_t, reg_lim, 2)
                        reg1_t, _ = mask_region(reg1_t, reg_lim, 2)
                        # print(reg1_t)
                        if len(reg1_t) > 0:
                            if type(reg1_t[0]) == list:
                                # reg1_t == [[num, num], [num, num]] or [[], [num, num]], etc

                                if type(reg1_t[1]) != list:
                                    print('error: list')

                                r = reg1_t
                                if r[0] != []:
                                    reg1_total[k] = r[0]
                                    reg1_t = r[0]
                                else:
                                    if debug:
                                        print('delete (1)')
                                    del reg1_total[k]
                                    # if len(reg1_total) == 0:
                                    #   break
                                    k -= 1
                                    break

                                if r[1] != []:
                                    reg1_total.append(r[1])
                                else:
                                    pass

                                # reg1_total[k] = reg1_t[0]
                                # reg1_total.append(reg1_t[1])
                                # reg1_t = reg1_t[0]

                            else:
                                # reg1_t == [num, num]
                                reg1_total[k] = reg1_t
                        else:
                            # reg1_t == []
                            if debug:
                                print('delete (2)')
                            del reg1_total[k]
                            # if len(reg1_total) == 0:
                            #   break
                            k -= 1
                            break
                k += 1

        else:
            # no limit
            if debug:
                print('no limit')


        if len(reg1_total) > 0:
            namelist[region_ids[i]] = reg1_total

        if debug:
            print(reg1, reg1_total)


    # debug checks
    regions_all = []
    for f in namelist:
        regions = namelist[f]
        for reg in regions:
            regions_all.append(reg)
    regions_all.sort()

    regions_all_cmb = combine_region(regions_all)

    if regions_all_cmb != regions_check:
        # combined regions start == combined regions end
        print('warning: regions_all != regions_check')

    for i in range(1,len(regions_all)):
        reg1 = regions_all[i-1]
        reg2 = regions_all[i]
        if reg1[1] != reg2[0]:
            b = 1
            # for reg_cmb in regions_all_cmb:
            for j in range(1,len(regions_all_cmb)):
                # if reg1[1] == reg_cmb[1] and reg2[0] == reg_cmb[0]:
                reg_cmb1 = regions_all_cmb[j-1]
                reg_cmb2 = regions_all_cmb[j]
                if reg1[1] == reg_cmb1[1] and reg2[0] == reg_cmb2[0]:
                    b = 0
            if b:
                """
                Final regions should make up a continuous range of time.
                If not, then the separation must be due to separations
                shown in the entire combined set, i.e.
                    region_ids = ['A','B','C']
                    regions    = [[10,20],[15,17],[30,40]]
                    P          = [1,2,3]
                    output:
                        namelist =
                        {'A': [[10,15], [17,20]],
                         'B': [[15,17]],
                         'C': [[30,40]]}
                        regions_all = 
                            [[10,15], [15,17], [17,20], [30,40]]
                        regions_all_cmb = 
                            [[10,20],[30,40]]
                    
                We should expect 10-15, 15-17, 17-20 to be continuous,
                then a separation only at the same place that there
                is separation in the combined set (regions_all_cmb),
                which is evident at 20-30.

                Floating/roundoff error may trigger this warning.
                """
                print(reg1, reg2)
                print('warning: reg')


    return namelist, regions_all


def find_region(t, dt, dt_buffer=0.0):

    """
    This function creates regions based on point-data. So,
    given an array t, each point will be assigned a region
    defined by -dt/+dt. These small regions are all collapsed
    into larger regions that define full 1d spaces of point
    data.

    Example:
        import numpy as np
        import icesatPlot as ip
        import region_detect as rd
        x = np.linspace(0,100)
        dx = rd.est_dx(x)
        x_cut, _ = rd.filt_reg_x(x, [[10,20], [30,40], [50,60]])
        regions, index = rd.find_region(x_cut, 2*dx)

        fig, ax = ip.make_fig()
        ax.plot(x_cut, x_cut, '.')
        ylim = ax.get_ylim()
        for reg in regions:
            ax.plot([reg[0], reg[0]], ylim, 'g--') # start
            ax.plot([reg[1], reg[1]], ylim, 'r--') # end
        fig.show()

    """

    degrade = []
    degrade_index_red = []
    if len(t) > 0:

        # This works by first giving every outlier point a region defined by
        # deg_definition_dt. Then, regions are added to degrade_index, and 
        # each overlapping super region is separated from one-another.
        # Then the first and last points of each super region are used to define
        # the overall degrade at that set of points.

        deg_regions = [[t[0] - dt, t[0] + dt]]
        degrade_index = [[0]]
        for k in range(1, len(t)):
            deg_t = t[k]
            deg_regions.append([deg_t - dt, deg_t + dt])

            # deg_start1 = deg_regions[k-1][0]
            deg_end1 = deg_regions[k-1][1]
            deg_start2 = deg_regions[k][0]
            # deg_end2 = deg_regions[k][1]

            if deg_end1 > deg_start2:
                # degrade regions overlap
                degrade_index[-1].append(k)
            else:
                degrade_index.append([])
                degrade_index[-1].append(k)


        # degrade_index[d] is a list of k-indices, the degrade regions
        # that all overlap, or a single region if len(degrade_index[d]) == 1

        for d in range(len(degrade_index)):
            # handles single point or multiple regions that overlap
            # in single point case, degrade_index[d][0] == degrade_index[d][-1]
            d_start = deg_regions[degrade_index[d][0]][0]   # "deg_start1", or "deg_start1" in single-point case
            d_end = deg_regions[degrade_index[d][-1]][1]    # "deg_end2", or "deg_end1" in single-point case
            d_start += (dt - dt_buffer)
            d_end += (-dt + dt_buffer)

            degrade.append([d_start, d_end])
            degrade_index_red.append([degrade_index[d][0], degrade_index[d][-1]])
            # fp_outlier_deg_h.write('%d %15.6f %15.6f %15.6f\n' % (i+1, degrade_start, degrade_end, degrade_end - degrade_start))

    return degrade, degrade_index_red


def get_total_region(deg1, deg2):
    # This function compares the degrade ranges and
    # returns a single total region
    
    # assumes start1 < end1, start2 < end2
    if deg1[1] < deg2[0]:
        # non-overlapping regions
        overlap = 0
        return deg1, overlap
    else:
        overlap = 1
        deg = np.array([deg1, deg2])
        c_reg = deg.T
        index_start = np.argmin(c_reg[0])
        index_end = np.argmax(c_reg[1])

        dr = [c_reg[2][index_start], c_reg[3][index_end]]

        deg_new = [c_reg[0][index_start], c_reg[1][index_end],
                   c_reg[2][index_start], c_reg[3][index_end]]

        return deg_new, overlap


def combine_region(degrades):

    """
    Combines a list of regions into a single region, if
    they overlap.

    degrades is a list of [start, end] degrades
        [[start, end]
         [start, end]
         [...]]

    Addition 9/24/19
    - degrades can also be a list of [start, end, dr1, dr2] regions
    - this enables each smaller region to encompass a larger zone,
      virtually, and the final combined region will only use the end-point
      dr's
    Format:
        [[start+dr1, end+dr2, dr1, dr2],
         [start+dr1, end+dr2, dr1, dr2],
         [...]]
    Example:
        [[0.0, 100.0, -10.0, 10.0],
         [-1.0, 50.0, -15.0, 25.0]]
        This would produce the combined region
        [[-1.0-15.0, 100.0+10.0]] == [[-16.0, 110.0]]

        Conceptually, if you had several regions that were separated by
        small gaps, you could combine them by artificially adding dr1
        and dr2.

    """

    # This works by choosing a first degrade,
    # then comparing the ones that come after it to see
    # if they're within the same degrade region.
    # Overlapping regions are combined into one degrade,
    # then a new "first" degrade is chosen, and so on.

    if type(degrades) != list:
        degrades = np.ndarray.tolist(degrades)
    
    n_deg = len(degrades)
    deg_reduced = [] # degrades with no overlapping regions
    if n_deg > 1:

        degrades_tr = np.transpose(degrades)
        if type(degrades_tr) == np.ndarray:
            if degrades_tr.dtype != 'O':
                length = len(degrades_tr)
                if length == 2:
                    degrades = [[deg[0], deg[1], 0, 0] for deg in degrades]
                elif length == 4:
                    pass
                else:
                    print('error: length of regions must be either 2 or 4')
                    return degrades
            else:
                print('error: inconsistent lengths in regions')
                return degrades

        degrades.sort() # sort by start time
        i = 0
        while i < n_deg:
            deg_new = degrades[i]
            j = i+1
            while j < n_deg:
                # recursively update the combined degrade
                deg_new, overlap = get_total_region(deg_new, degrades[j])
                if not overlap:
                    break
                j += 1

            # append the final modification
            deg_reduced.append(deg_new)
            i = j

        if length == 2:
            deg_reduced = [[deg[0], deg[1]] for deg in deg_reduced]

    elif n_deg == 1:
        deg_reduced = [degrades[0]]

    return deg_reduced