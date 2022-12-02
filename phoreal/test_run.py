import core


prDS = core.PhorealDataSet()

working_dir = '/home/mjh5468/ContainerDevelopment/data/icesat2_granules/'

prDS.load_directory(directory=working_dir)

prDS.build_granules_metadata()


print(prDS.granule_list)