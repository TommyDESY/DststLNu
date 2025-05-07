"""
Steering file for inclusive B->Xlnu decays at generator level.
Needed to generate the hybrid weights and estimate the efficiency.

Should only be run on the B->Xulnu signal MC samples WITHOUT SKIM!
"""

import basf2
import modularAnalysis as ma
from variables import variables as vm

from numpy import mask_indices
from ROOT import Belle2

import yaml

import argparse
import glob
import pathlib
from typing import List, Optional


def parse_cmdline_args() -> argparse.Namespace:
    """
    Creates the argument parser and returns command line arguments.
    """
    parser = argparse.ArgumentParser(description='Inclusive B->Xlv reconstruction steering file.')
    parser.add_argument(
        '-input',
        '--input-path',
        dest='input_path',
        help='Path to the file that has the input files stored',
    )

    parser.add_argument(
        '-output',
        '--output-file',
        dest='output_file',
        help='Full path of the output file',
    )
    return parser.parse_args()


def main():
    """
    This should be seen as a test.
    The correct way to process the sample is via the `steering_pipeline.py`.
    """

    args = parse_cmdline_args()
    # Input files
    input_files = [args.input_path]
    # Output file
    output_filename = f'{args.output_file}' if args.output_file else 'b2xlv_output_generator.root'
    # Create output path dirs
    output_split = output_filename.split('/')
    output_dir = '/'.join(output_split[:-1])
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    # User defined variables
    basf2.add_module_search_path('.')
    path = basf2.create_path()
    basf2.register_module('InclusiveVariables')  # now you can use it
    path.add_module('InclusiveVariables', shared_lib_path='./libbasf2_patches.so')

    path.add_path(create_path(
        input_files,
        output_path=output_filename,
    ))
    basf2.process(path, 100000)
    print(basf2.statistics)
    return


def get_gen_variables():
    """
        gets the list of generator level variables. (Subset of inclusive_variables.py (alternatively somehow flag them there.))
    """
    upsilon_variables = [
        'eventRandom',
        ('gen_lep_p_B', 'genPstarL'),
        ('gen_lep_E_B', 'genEstarL'),
        'genPplus',
        ('genq2', 'genqSquared'),
        'genMx',
        'genEx',
        'genPXx',
        'genPYx',
        'genPZx',
        'genMx2',
        'genW',
        'genCosThetaL',
        'genCosThetaV',
        'genChiPrime',
        'genTagBPx',
        'genTagBPy',
        'genTagBPz',
        'genTagBE',
        ('Btag_gen_PDG', 'daughter(0, mcPDG)'),
        ('Btag_mdstIndex', 'daughter(0, mdstIndex)'),
        (
            'genY4S_nPi0',
            'genUpsilon4S(nMCDescendentsWithPDGCode(111))',
        ),
        (
            'genY4S_nPion',
            'genUpsilon4S(nMCDescendentsWithPDGCode(211))',
        ),
        (
            'genY4S_nKaon',
            'genUpsilon4S(nMCDescendentsWithPDGCode(321))',
        ),
        (
            'genBsig_nPi0',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(111)))',
        ),
        (
            'genBsig_nPiCh',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(211)))',
        ),
        (
            'genBsig_nKaonL',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(130)))',
        ),
        (
            'genBsig_nKaonS',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(310)))',
        ),
        (
            'genBsig_nKaon0',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(311)))',
        ),
        (
            'genBsig_nKaonCh',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(321)))',
        ),
    ]

    lepton_variables = [
        ('lep_mdstIndex', 'daughter(1, mdstIndex)'),
        ('lep_motherMdstIndex', 'daughter(1, mcMother(mdstIndex))'),
        ('lep_PDG', 'daughter(1, PDG)'),
        ('lep_mcPDG', 'daughter(1, mcPDG)'),
        ('Bsig_gen_PDG', 'daughter(1, mcMother(PDG))'),
        ('Y4S_gen_PDG', 'daughter(1, mcMother(mcMother(PDG)))'),
        (
            'Upsilon4S_mcMother_gen_PDG',
            'daughter(1, mcMother(mcMother(mcMother(PDG))))',
        ),
        # these are for matching the generated decay process
        ('X_gen_PDG', 'daughter(1, mcMother(mcDaughter(0, PDG)))'),
        (
            'Xdecay0_gen_PDG',
            'daughter(1, mcMother(mcDaughter(0, mcDaughter(0, PDG))))',
        ),
        (
            'Xdecay1_gen_PDG',
            'daughter(1, mcMother(mcDaughter(0, mcDaughter(1, PDG))))',
        ),
        (
            'Xdecay2_gen_PDG',
            'daughter(1, mcMother(mcDaughter(0, mcDaughter(2, PDG))))',
        ),
        (
            'Xdecay3_gen_PDG',
            'daughter(1, mcMother(mcDaughter(0, mcDaughter(3, PDG))))',
        ),
        (
            'Xdecay4_gen_PDG',
            'daughter(1, mcMother(mcDaughter(0, mcDaughter(4, PDG))))',
        ),
        (
            'XGD0_gen_PDG',
            'daughter(1, mcMotherSkipW(mcDaughterSkipW(0, mcDaughterSkipW(0, mcDaughterSkipW(0, PDG)))))',
        ),  # X granddaughter 1 gen PDG
        (
            'XGD1_gen_PDG',
            'daughter(1, mcMotherSkipW(mcDaughterSkipW(0, mcDaughterSkipW(0, mcDaughterSkipW(1, PDG)))))',
        ),  # X granddaughter 2 gen PDG
        ('lep_gen_PDG', 'daughter(1, mcMother(mcDaughter(1, PDG)))'),
        ('neutrino_gen_PDG', 'daughter(1, mcMother(mcDaughter(2, PDG)))'),
        ('Bdecay3_gen_PDG', 'daughter(1, mcMother(mcDaughter(3, PDG)))'),
        ('Bdecay4_gen_PDG', 'daughter(1, mcMother(mcDaughter(4, PDG)))'),
        ('Bdecay5_gen_PDG', 'daughter(1, mcMother(mcDaughter(5, PDG)))'),
        ('Bdecay6_gen_PDG', 'daughter(1, mcMother(mcDaughter(6, PDG)))'),

        # generator level quantities. Operate on the lepton to get the correct MC particles.
        'genLepBFrameE',
        'genLepBFramePx',
        'genLepBFramePy',
        'genLepBFramePz',
        'genNeuBFrameE',
        'genNeuBFramePx',
        'genNeuBFramePy',
        'genNeuBFramePz',
        'genXBFrameE',
        'genXBFramePx',
        'genXBFramePy',
        'genXBFramePz',

        # Hadron daughter generator level information. Needed for Hammer FF reweighting.
        'genXD1BFrameE',
        'genXD1BFramePx',
        'genXD1BFramePy',
        'genXD1BFramePz',
        'genXD2BFrameE',
        'genXD2BFramePx',
        'genXD2BFramePy',
        'genXD2BFramePz',
        'genXD3BFrameE',
        'genXD3BFramePx',
        'genXD3BFramePy',
        'genXD3BFramePz',
        'genXGD1BFrameE',
        'genXGD1BFramePx',
        'genXGD1BFramePy',
        'genXGD1BFramePz',
        'genXGD2BFrameE',
        'genXGD2BFramePx',
        'genXGD2BFramePy',
        'genXGD2BFramePz',

        # Pion and Kaon multiplicities on the signal side
        (
            'genBsig_nPiZero',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(111)))',
        ),
        (
            'genBsig_nPiCh',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(211)))',
        ),
        (
            'genBsig_nKaonCh',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(321)))',
        ),
        (
            'genBsig_nKaonS',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(310)))',
        ),
        (
            'genBsig_nKaonL',
            'daughter(1, mcMother(nMCDescendentsWithPDGCode(130)))',
        ),
    ]

    var_list = []
    for var in [*lepton_variables, *upsilon_variables]:
        if type(var) == tuple:  # change this
            alias, name = var
            vm.addAlias(alias, name)
            var_list.append(alias)
        else:
            var_list.append(var)
    return var_list


def create_path(
    input_files: List[str],
    output_path: Optional[str] = None,
) -> 'basf2.Path':
    """
    Creates the `basf2` path object to be executed for the analysis.

    :return: Analysis `basf2` path object.
    :rtype: basf2.Path
    """

    main_path = basf2.create_path()

    ma.inputMdstList(
        filelist=input_files,
        entrySequences=[],
        path=main_path,
    )

    # load particles
    ma.fillParticleListFromMC('B+:gen', '', addDaughters=True, path=main_path)
    ma.fillParticleListFromMC('B0:gen', '', addDaughters=True, path=main_path)

    ma.fillParticleListFromMC('mu+:gen',
                              'abs(mcMother(PDG))==511 or abs(mcMother(PDG))==521',
                              addDaughters=True,
                              path=main_path)
    ma.fillParticleListFromMC('e+:gen',
                              'abs(mcMother(PDG))==511 or abs(mcMother(PDG))==521',
                              addDaughters=True,
                              path=main_path)

    common_cut = ''
    ma.reconstructDecay('Upsilon(4S):gen_charged_mu -> B+:gen mu-:gen',
                        common_cut,
                        allowChargeViolation=True,
                        path=main_path)
    ma.reconstructDecay('Upsilon(4S):gen_charged_e -> B+:gen e-:gen',
                        common_cut,
                        allowChargeViolation=True,
                        path=main_path)

    ma.reconstructDecay('Upsilon(4S):gen_mixed_mu -> B0:gen mu-:gen',
                        common_cut,
                        allowChargeViolation=True,
                        path=main_path)
    ma.reconstructDecay('Upsilon(4S):gen_mixed_mu_ws -> B0:gen mu+:gen',
                        common_cut,
                        allowChargeViolation=True,
                        path=main_path)
    ma.reconstructDecay('Upsilon(4S):gen_mixed_e -> B0:gen e-:gen',
                        common_cut,
                        allowChargeViolation=True,
                        path=main_path)
    ma.reconstructDecay('Upsilon(4S):gen_mixed_e_ws -> B0:gen e+:gen',
                        common_cut,
                        allowChargeViolation=True,
                        path=main_path)

    ma.copyLists('Upsilon(4S):gen_combined', [
        'Upsilon(4S):gen_charged_mu',
        'Upsilon(4S):gen_charged_e',
        'Upsilon(4S):gen_mixed_mu',
        'Upsilon(4S):gen_mixed_mu_ws',
        'Upsilon(4S):gen_mixed_e',
        'Upsilon(4S):gen_mixed_e_ws',
    ],
                 path=main_path)

    ### Variables and Ntuple creation
    variables = get_gen_variables()
    ma.variablesToNtuple(
        'Upsilon(4S):gen_combined',
        variables=variables,
        filename=output_path,
        treename='InclusiveSLGen',
        path=main_path,
    )
    return main_path


if __name__ == '__main__':
    main()
