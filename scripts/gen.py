import basf2 as b2
import generators as ge
import modularAnalysis as ma
import sys

def main(argv):

    inputfile_name = sys.argv[1]
    my_path = b2.create_path()

    # Setting up number of events to generate
    ma.setupEventInfo(noEvents=50000, path=my_path)

    # Adding generator
    ge.add_evtgen_generator(
        path=my_path,
        finalstate='signal',
        signaldecfile=b2.find_file(f'{inputfile_name}.dec'),
        )

    # If the simulation and reconstruction is not performed in the sam job,
    # then the Gearbox needs to be loaded with the loadGearbox() function.
    #ma.loadGearbox(path=my_path)

    my_path.add_module(
        'RootOutput', 
        outputFileName= f'root/{inputfile_name}.root',
        )

    b2.process(my_path)

if __name__ == "__main__":
    main(sys.argv[1:])
