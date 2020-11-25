

import os, sys
import numpy as np
import icesatUtils as iu
import h5py as h5


def search(DIR, search_list, logical='and', case=True, ignore_perms=True, ext=None, recursive=True, debug=0):

    """
    Finds files in DIR, that have the criteria in search_list. The logical specifies whether search_list 
    entries are 'or' or 'and'. case turns on case-sensitivity. ignore_perms... ignores permissions. ext
    is there if you'd like to look for a specific kind of file (typically h5, but it could be anything).
    recursive means it'll search through all sub-directories, so it could take some time when this is True.

    Input:
        DIR - top directory to search recursively through (if recursive=True)
        search_list - substrings to search for, such as ['ATL03', 'ATL08'] or even ['20181016000635']
        logical - 'or' or 'and' for pairing search_list entries
        case - in case upper/lowercase is necessary
        ignore_perms - .. ignore permissions catches, which stops the code
        ext - file extension filter, so if ext='h5' it'll only find h5 files (case-sensitive regardless of case bool)
        recursive - if true, code searches all sub-directories (default)

    Output:
        list of string full-path files

    Example:
        # find all files in bigtex data that contain ATL03 -or- ATL08 (default is 'and')
        import file_search as fs
        files = fs.search(DATA_DIR, ['ATL03', 'ATL08'], logical='or')

    """

    if debug:
        if ignore_perms:
            print('ignoring permissions denied files')
        else:
            print('checking permissions denied files')


    logical = logical.lower()
    if not (logical == 'or' or logical == 'and'):
        print('error: logical argument = %s; must be either \'or\' or \'and\'' % logical)
        return []
        
    if len(search_list) == 0:
        return []

    import subprocess
    cmd = 'find %s' % DIR
    if not recursive:
        cmd += ' -maxdepth 1'
    for i, s in enumerate(search_list):
        if not case:
            s = s.upper()
        if i == 0:
            cmd += ' -name \'*%s*\' -type f' % s
        else:
            cmd += ' -%s -name \'*%s*\' -type f' % (logical, s)

    # if ignore_perms:
    cmd += " 2>&1 | grep -v \"Permission denied\"" # ignore permission denied files
    if debug:
        print(cmd)

    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        return []

    s = output.decode()
    files = s.split('\n')[:-1]

    if ext != None:
        files = filter_ext(files, ext)

    return files



def show_h5(h5_file, key_main='/'):

    """
    Input h5 file, outputs groups/datasets
    """
    fp = h5.File(h5_file, 'r')

    # key_main = '/gt2r'
    for key in fp[key_main].keys():
        print(key_main + '/' + key)
        try:
            for val in list(fp[key_main + '/' + key].keys()):
                print(val)
        except:
            pass
        print('\n')

    fp.close()

        
def filter_ext(files, ext='h5'):
    """
    This filter returns only .h5 files (by default), given any kind of files
    otherwise, set ext to the desired extension to filter out
    Input:
        list of string full-path files, extension to filter for

    Output:
        list of string filtered full-path files

    Example:
    # filter for h5 files only
    import file_search as fs
    files_filt = fs.filter_ext(files, ext='h5')

    """
    # files_ext = list(filter(lambda x: x.split('/')[-1][-len(ext):] == ext, files))
    files_ext = list(filter(lambda x: os.path.basename(x)[-len(ext):] == ext, files))
    return files_ext



def filter_atl(files, year=None, doy=None, hms=None, track=None, cycle=None, f_type=None,
                     release=None, version=None, name=None, debug=0):

    """
    Input:
        All inputs can be stacked, such as filter_atl(files, year=2020, doy=250, cycle=506)

        files - list of full-path atl files, mixed 03, 08, etc is ok.
                Could use non-ATL files, but it will likely return []
                This defaults to 'and' logic, so to get 'or' (such as
                ATL files with this year -or- ATL files with this 
                track), you'll have to run the function twice with
                two different filters (in this case, year and track).

        year - year, either string or int
        doy - day of year, either string or int
        hms - string of hour, minute, sec of atl file
        track - track number, either string or int
        cycle - cycle number, either string or int
        f_type - either 'rapid' or 'final'
        release - release number, string or int
        version - version number, string or int
        name - any string to match part of the file to,
                i.e. filter_atl(files, name='202002') finds
                        files with year 2020 and month 02
        debug - in case user wants warning messages

    Output:
        filtered list of files

    Example:
        import file_search as fs
        DIR = fs.get_dir('003', 9)

        # get all release 3 atl09 files
        files_09 = fs.search(DIR, ['ATL09'], ext='h5', recursive=False)

        # filter for year 2019, month 8, day 15, cycle 401
        files_filt = fs.filter_atl(files_09, name='20190815', cycle=401)

           
    If you want to add any args to this, you'll have to modify
        args
        args_desc
        get_h5_meta()
        maybe apply_filter() in this function, if get_h5_meta
            needs special adapation (see year/doy)
    """

    if year != None:
        year = str(year)

    if doy != None:
        doy = str(doy).zfill(3)

    if hms != None:
        hms = str(hms)

    if track != None:
        track = int(track)

    if release != None:
        # release = str(release).upper()
        release = str(release).zfill(3)

    if version != None:
        version = str(version).zfill(2)

    if cycle != None:
        cycle = int(cycle)

    if f_type != None:
        f_type = f_type.lower()

    if name != None:
        name = str(name)


    args = [year, doy, hms, track, cycle, f_type, release, version, name]
    args_desc = ['year', 'doy', 'hms', 'track', 'cycle', 'f_type', 'release', 'version', 'name']
    


    if len(files) == 0:
        if debug:
            print('warning: len(files) == 0)')
        return files

    if (np.array(args) == None).all():
        if debug:
            print('warning: no filter entered')
        return files

    else:

        filter_dict = {}
        for i, arg in enumerate(args):
            if arg != None:
                filter_dict[args_desc[i]] = arg


        def apply_filter(files, f_desc, f_val):

            def warn_msg(f):
                print('warning: could not get meta from %s' % f)

            files_new = []
            if f_desc == 'year':
                for file_full in files:
                    try:
                        y, d = iu.get_h5_meta(file_full, rtn_doy=True, file_start='ATL')
                    except ValueError:
                        if debug:
                            warn_msg(file_full)
                        continue
                    if y == f_val:
                        files_new.append(file_full)
                return files_new

            elif f_desc == 'doy':
                for file_full in files:
                    try:
                        y, d = iu.get_h5_meta(file_full, rtn_doy=True, file_start='ATL')
                    except ValueError:
                        if debug:
                            warn_msg(file_full)
                        continue
                    if d == f_val:
                        files_new.append(file_full)
                return files_new

            elif f_desc == 'hms':
                for file_full in files:
                    try:
                        h,m,s = iu.get_h5_meta(file_full, f_desc, file_start='ATL')
                        hms_rtn = h+m+s
                    except ValueError:
                        if debug:
                            warn_msg(file_full)
                        continue
                    if hms_rtn == f_val:
                        files_new.append(file_full)
                return files_new

            elif f_desc == 'f_type':
                for file_full in files:
                    try:
                        f_type_rtn = iu.get_h5_meta(file_full, f_desc, file_start='ATL')
                    except ValueError:
                        if debug:
                            warn_msg(file_full)
                        continue
                    if f_type_rtn[0] == f_val[0]:
                        files_new.append(file_full)
                return files_new

            elif f_desc == 'name':
                for file_full in files:
                    # file = file_full.split('/')[-1]
                    file = os.path.basename(file_full)
                    if name in file:
                        files_new.append(file_full)
                return files_new

            else:
                # track, cycle, release, version
                for file_full in files:
                    try:
                        val = iu.get_h5_meta(file_full, f_desc, file_start='ATL')
                    except ValueError:
                        if debug:
                            warn_msg(file_full)
                    if val == f_val:
                        files_new.append(file_full)
                return files_new


        for f_desc in filter_dict:
            f_val = filter_dict[f_desc]
            files = apply_filter(files, f_desc, f_val)


        return files


# def check_file_start(file_sub, file_start, debug=0):
#     i0 = 0
#     try:
#         i0 = file_sub.index(file_start)
#     except ValueError:
#         if debug:
#             print('warning: substring %s not found in %s' % (file_start, file_sub))
#     return i0

def remove_duplicate_atl(files_h5, release='latest', version='latest', file_start='ATL', debug=0):

    """
    This function takes in a list of files, and removes
    files that are very similar by prioritizing release,
    then version. Any duplication at the end should be
    based purely on differences in _sub or _nsidc naming
    convention.

    This function does not prioritize by file size - this
    would be a good feature to add.

    Input:
        files_h5 - a list of files, typically similar but not necessarily
        release - defined as 'latest', but could be any release string 
                    such as '001', '002', etc
                    Code will not throw a warning if release isn't found,
                    it'll just return []
        version - similar to release, either 'latest' or '01', '02', etc
        file_start - if error_processing prefix is on the ATL file
        debug - for user warnings

    Output:
        return files_h5_new
        A list of files with redundant ones removed

    Example:
        see if __name__ == 'main': code
        and run
        $ python file_system.py

    """
    def check_file_start(file_sub, file_start, debug=0):
        i0 = 0
        try:
            i0 = file_sub.index(file_start)
        except ValueError:
            if debug:
                print('warning: substring %s not found in %s' % (file_start, file_sub))
        return i0

    if release != 'latest':
        release = str(release).zfill(3)
    if version != 'latest':
        version = str(version).zfill(2)

    # this removes exact duplicate names
    # files_h5_sub = [file_h5.split('/')[-1][:36] for file_h5 in files_h5]
    files_h5_sub = []
    for file_h5 in files_h5:
        # file_h5_sub = file_h5.split('/')[-1]
        file_h5_sub = os.path.basename(file_h5)
        i0 = check_file_start(file_h5_sub, file_start, debug)
        file_h5_sub_match = file_h5_sub[i0:i0+36]
        files_h5_sub.append(file_h5_sub_match)

    _, index = np.unique(files_h5_sub, return_index=True)
    files_h5 = list(np.array(files_h5)[index])

    # this next bit prioritizes either release or version, if applicable
    files_h5_new = []
    for file0 in files_h5:
        # files_new = filter_atl(files_h5, name=file0)

        # file0_sub = file0.split('/')[-1]
        file0_sub = os.path.basename(file0)
        i0 = check_file_start(file0_sub, file_start, debug)
        files_name = filter_atl(files_h5, name=file0_sub[i0:i0+29])
        # print(i0)
        # print(files_name)
        if len(files_name) > 1:
            # multiple versions or releases, or some 
            # files have nreq or sub or some other bs

            if release == 'latest':
                # filter by release first
                r = []
                for file_n in files_name:
                    # file_n_sub = file_n.split('/')[-1]
                    file_n_sub = os.path.basename(file_n)
                    i0 = check_file_start(file_n_sub, file_start, debug)
                    r.append(int(file_n_sub[i0+30:i0+33]))
                r_max = max(r)
                files_name_r = filter_atl(files_name, release=r_max)

            else:
                files_name_r = filter_atl(files_name, release=release)

            if version == 'latest':
                # filter by version second
                v = []
                for file_n in files_name_r:
                    # file_n_sub = file_n.split('/')[-1]
                    file_n_sub = os.path.basename(file_n)
                    i0 = check_file_start(file_n_sub, file_start, debug)
                    v.append(int(file_n_sub[i0+34:i0+36]))
                v_max = max(v)
                files_name_v = filter_atl(files_name_r, version=v_max)

            else:
                files_name_v = filter_atl(files_name_r, version=version)

            files_h5_new.append(files_name_v[0])

        else:
            files_h5_new.append(files_name[0])

    files_h5_new = list(np.unique(files_h5_new))
    return files_h5_new


def find_match(file, search_type, SEARCH_DIR=None, recursive=True, match_ymd=True, match_trackcycle=True,
                 match_release=False, match_version=False, debug=0):

    """
    This code finds matches between ATL files, 03 to 08, 03 to 09, 09 to 12 - anything.

    Input:
        file - full-path ATL file (such as 03)
        search_type - must be of format "ATL##", not case-sensitive
        SEARCH_DIR - optional, speeds up processing if chosen
        recursive - sets search to look thru dirs recursively or not
        match_ymd - match year, month, day (default)
        match_trackcycle - match both track and cycle (default)
        match_release - match release number
        match_version - match version number

    Output:
        return files_out
        a list of search_type "ATL##" files that match the given file

    Caveats:
        NOTE
        The boolean matches occur via priority, so you could not set
        match_version = True, then set match_ymd = False; it would just
        match all the way -up to- the version number, including y,m,d.

        This code may find several files that correspond to a single input
        file to match. For this, see remove_duplicate_atl().

    Example:
        import file_search as fs
        file_03 = '/bigtex_data/data/release/002/ATL03_r002/Finland/ATL03_20190819234642_08080403_002_01_sreq_3040.h5'
        files_08 = fs.find_match(file_03, search_type='ATL08', SEARCH_DIR='/bigtex_data/data/release/002/ATL08_r002')

        file_08 = '/bigtex_data/data/release/002/ATL08_r002/ATL08_20190819234642_08080403_002_02.h5'
        files_03 = fs.find_match(file_08, search_type='ATL03', SEARCH_DIR='/bigtex_data/data/release/002/ATL03_r002')

    """

    file_start='ATL'
    search_type = search_type.lower()

    if SEARCH_DIR == None:
        SEARCH_DIR = get_dir('data')
        if debug:
            print('searching all of %s...'% SEARCH_DIR)


    # file_split = file.split('/')
    # file_sub = file_split[-1]
    file_sub = os.path.basename(file)
    index = file_sub.index(file_start)
    file_sub = file_sub[index:]
    
    # if match_release and match_version:
    if match_version:
        file_cmp = file_sub[:36] # version

    elif match_release:
        file_cmp = file_sub[:33] # release
    
    elif match_trackcycle:
        file_cmp = file_sub[:29] # track/cycle

    elif match_ymd:
        file_cmp = file_sub[:20] # ymd

    else:
        file_cmp = file_sub[:14] # hms

    file_cmp = list(file_cmp)

    file_cmp_new = file_cmp
    search_type_numeral = search_type[3:5]
    file_cmp_new[3] = search_type_numeral[0]
    file_cmp_new[4] = search_type_numeral[1]
    file_cmp_new = ''.join(file_cmp_new)

    if debug:
        print(SEARCH_DIR)
        print(file_cmp_new)

    files_out = search(SEARCH_DIR, [file_cmp_new], ext='h5', recursive=recursive)

    if len(files_out) == 0:
        if debug:
            print('no files found')

    return files_out


def get_dir(DIR='data', atl=3): #, year=None, doy=None):

    """
    Exclusively gets data directory (relative to user) for
    bigtex.
    Example:
        # release 003 for atl03
        import file_search as fs
        DATA_DIR = fs.get_dir('003', 3)
        files_03 = fs.search(DATA_DIR, ['ATL03'], ext='h5')
        
    """

    BIGTEX_DIR = iu.get_root_dir('data')

    if DIR == '':
        print('error: DIR must be defined for bigtex')
        return None

    DIR = DIR.lower()
    if DIR == 'data':
        # return BIGTEX_DIR + '/data/release'
        return os.path.join(BIGTEX_DIR, 'data', 'release')

    else:
        # return BIGTEX_DIR + '/data/release/%s/ATL%s_r%s' % (DIR.upper(), str(atl).zfill(2), DIR.upper())
        rel, d = DIR.upper(), 'ATL%s_r%s' % (str(atl).zfill(2), DIR.upper())
        return os.path.join(BIGTEX_DIR, 'data', 'release', rel, d)


if __name__ == "__main__":

    """
    This example matches all 03 files to all 08
    files in release 003
    """
    print('find atl03 and atl08 files for release 003')
    DIR_03 = get_dir('003', 3)
    DIR_08 = get_dir('003', 8)
    files_03 = search(DIR_03, ['ATL03'], ext='h5', recursive=False)
    # files_08 = search(DIR_08, ['ATL08'], ext='h5', recursive=False)
    
    print('match 03 files to 08 files')
    files_03_match = []
    files_08_match = []
    n_03 = len(files_03)

    from tqdm import tqdm
    for i in tqdm(range(n_03)):
        file_03 = files_03[i]
        files_08 = find_match(file_03, 'atl08', DIR_08) #, recursive=False)
        # files_09 = find_match(file_03, 'atl09', DIR_09) #, recursive=False)
        if len(files_08) == 0:
            continue
        if len(files_08) > 1:
            files_08 = remove_duplicate_atl(files_08)
            # for file_08 in files_08:
            #     print(file_08)
            # print('')
            # iu.pause()

        file_08 = files_08[0]
        files_03_match.append(file_03)
        files_08_match.append(file_08)

    #     if len(files_03_match) > 10:
    #         break


    # if 1:
    #     import icesatReader as ir
    #     # file_03, file_08 = files_03_match[0], files_08_match[0]
    #     num_files = len(files_03_match)
    #     for f in range(num_files):
    #         print(f)
    #         file_03 = files_03_match[f]
    #         file_08 = files_08_match[f]
    #         try:
    #             atl03 = ir.get_atl03_struct(file_03, 'gt1r', file_08)
    #         except:
    #             continue

    #         if not atl03.dataIsMapped:
    #             continue

    #         ir.write_atl03_las(atl03, '/LIDAR/server/poseidon_files/USERS/jsipps/scripts_lidar/icesat_analysis')
    #         break