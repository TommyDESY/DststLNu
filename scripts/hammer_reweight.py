import numpy as np
import pandas as pd
from typing import Dict, Optional, Union
import copy
import logging
from hammer.hammerlib import FourMomentum, Hammer, Particle, Process
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))
import database
from database.utility import mc_pdg_codes, pdg_to_string_dict


def add_ff_vars_to_hammer(
    HAMMER,
    model_dict: dict,
    dec_dict: dict,
    new_dict: Optional[dict] = None,
    verbose: bool = True,
):
    """
    Function that initializes the FFs of a specific decay in the Hammer class.
    :param HAMMER: Hammer class to be used.
    :param model_dict: Dictionary of all FF models that are used for all included decays.
    :param dec_dict: Dictionary of the FF model that was used in the dec file.
    :param new_dict: Dictionary of the FF model to be used now. If None, the dec_dict is used.
    :param verbose: If true, the different steps are printed.
    """
    if new_dict is None:
        new_dict = dec_dict

    allowed_decays = [f"B{dec_dict['Xc_string']}{dec_dict['lep_string']}Nu"]

    HAMMER.include_decay(allowed_decays)
    if verbose:
        logging.info(f'Added the following decays: {allowed_decays}')

    # No need to set options if we start from ISGW2
    # HAMMER.set_options(f"Bto{dec_dict['Xc_string']}{dec_dict['model']}: {dec_dict['param_initializer']}" % tuple([
    #     dec_dict[f'param_{i}'] if dec_dict[f'name_{i}'] != 'DelMbc' else dec_dict['param_mc']
    #     for i in range(dec_dict['num_params'])
    # ]))

    model_dict[f"B{new_dict['Xc_string']}"] = f"{new_dict['model']}_norm"
    HAMMER.add_ff_scheme(f"Bto{new_dict['Xc_string'].replace('*', 'st')}_{new_dict['model']}_norm", model_dict)
    HAMMER.set_options(f"Bto{new_dict['Xc_string']}{new_dict['model']}_norm : {new_dict['param_initializer']}" % tuple([
        new_dict[f'param_{i}'] if new_dict[f'name_{i}'] != 'DelMbc' else new_dict['param_mc']
        for i in range(new_dict['num_params'])
    ]))

    if verbose:
        logging.info(40 * '*')
        logging.info(f"Added decay B{dec_dict['Xc_string']}{dec_dict['lep_string']}Nu")
        if dec_dict['model'] == 'ISGW2':
            logging.info(
                f"Added the following FF input: Bto{dec_dict['Xc_string']}{dec_dict['model']}. No parameters are changed from default."
            )
        else:
            pass
            # logging.info(
            #     f"Added the following FF input: Bto{dec_dict['Xc_string']}{dec_dict['model']}: {dec_dict['param_initializer']}"
            #     % tuple([
            #         dec_dict[f'param_{i}'] if dec_dict[f'name_{i}'] != 'DelMbc' else dec_dict['param_mc']
            #         for i in range(dec_dict['num_params'])
            #     ]))

        logging.info(
            f"Added the new FF model: Bto{new_dict['Xc_string']}{new_dict['model']}_{dec_dict['lep_string']}norm: {new_dict['param_initializer']}"
            % tuple([
                new_dict[f'param_{i}'] if new_dict[f'name_{i}'] != 'DelMbc' else new_dict['param_mc']
                for i in range(new_dict['num_params'])
            ]))

    return HAMMER

def initialize_BtoXcLepNu_hammer(
    ffs: Dict[str, Union[int, float, np.array]],
    mode='XcEllNu',
    verbose: bool = True,
):
    """
    Function that implements the BtoXcEllNu and BtoXcTauNu form factor weights.
    :param HAMMER: Hammer class to be used. If None, a new class is created.
    :param mode: Mode to specify which models should be added. Choose b/w XcEllNu, XcTauNu or
                 individual decays. Note that the latter case will probably be the fastest.
                 XcEllNu and XcTauNu decays should not necesarily be treated together as the
                 B -> D(*/**) decays with an Ell would be affected by the model of this decay
                 together with a Tau as well.
    :param verbose: If true, this comes with printouts.
    """

    HAMMER = Hammer()

    model_dict = {}
    dicts = [
        ffs['BtoD0stEllNu_new'], ffs['BtoD0stEllNu_new'],
    ]
    for dct in dicts:
        model_dict[f"B{dct['Xc_string']}"] = dct['model']
        # in the narrow D** cases, I need to add the D** decay FF scheme as well (according to Dean):
        # otherwise, D**2*->D breaks. Moreover, it changes the outcome of D**2*/D**1->D*!
        if dct['Xc_string'] in ['D**1', 'D**2*']:
            # model_dict[f"{dct['Xc_string']}D*Pi"] = 'PW'
            # if dct['Xc_string'] == 'D**2*':
            model_dict[f"{dct['Xc_string']}DPi"] = 'PW'
    if verbose:
        logging.info('Model dict:')
        logging.info(model_dict)
    # tell hammer which form factor model was used in the MC (this can be found
    # in the decay files of Belle II):
    HAMMER.set_ff_input_scheme(model_dict)

    # Note: we can't use the BGLB2Dean weights anymore as it's become impossible to install Hammer on NAF
    # Instead we use the BGL model with specific input parameters that match the BGLB2Dean model
    # B0 -> D1* (-> D pi) l nu
    HAMMER = add_ff_vars_to_hammer(HAMMER,
                                   model_dict=copy.deepcopy(model_dict),
                                   dec_dict=ffs['BtoD0stEllNu_new'],
                                   new_dict=ffs['BtoD0stEllNu_new'],
                                   verbose=verbose)  # ISGW2 -> BLR

    HAMMER.set_units('GeV')  # this is also the default
    logging.info('units')
    HAMMER.init_run()
    logging.info('init_run')
    if verbose:
        logging.info(40 * '-')
        logging.info('HAMMER successfully initialized! The form factor schemes are:')
        logging.info([schemename for schemename in HAMMER.get_ff_scheme_names()])
        logging.info(40 * '-')
    return HAMMER


def proc_BtoXcLepNu_event(BPDG,
                          XcE,
                          Xcpx,
                          Xcpy,
                          Xcpz,
                          XcPDG,
                          XcD1E,
                          XcD1px,
                          XcD1py,
                          XcD1pz,
                          XcD1PDG,
                          XcD2E,
                          XcD2px,
                          XcD2py,
                          XcD2pz,
                          XcD2PDG,
                          XcD3E,
                          XcD3px,
                          XcD3py,
                          XcD3pz,
                          XcD3PDG,
                          XcGD1E,
                          XcGD1px,
                          XcGD1py,
                          XcGD1pz,
                          XcGD1PDG,
                          XcGD2E,
                          XcGD2px,
                          XcGD2py,
                          XcGD2pz,
                          XcGD2PDG,
                          lepE,
                          leppx,
                          leppy,
                          leppz,
                          lepPDG,
                          BD3E,
                          BD3px,
                          BD3py,
                          BD3pz,
                          BD3PDG,
                          HAMMER,
                          weightdict,
                          ratedict,
                          BE=1,
                          Bpx=0,
                          Bpy=0,
                          Bpz=0,
                          ignore_photons=False,
                          verbose=False):

    if verbose:
        logging.info('Starting...')
    if abs(BPDG) == 521:
        BE = 5.27934 # hardcoded from https://gitlab.desy.de/belle2/software/basf2/-/blob/main/framework/particledb/data/evt.pdl#L448
    elif abs(BPDG) == 511:
        BE = 5.27965 # hardcoded from https://gitlab.desy.de/belle2/software/basf2/-/blob/main/framework/particledb/data/evt.pdl#L450
    schemename_strings = np.array([schemename for schemename in HAMMER.get_ff_scheme_names()])
    #for schemename in schemename_strings:
    #    if schemename not in weightdict.keys():
    #        raise ValueError(f"ERROR! {schemename} is part of HAMMER but not of the weightdict...")
    HAMMER.init_event()  # in here, the process is reset that might have been stored previously in the class
    PROC = Process()  # a new process needs to be used for every event, otherwise I run into a secret infinity-loop!
    # add all particles that participate in the decay
    # B and D*:
    if verbose:
        logging.info(f'B: {BPDG}')
    BID = PROC.add_particle(Particle(FourMomentum(BE, Bpx, Bpy, Bpz), int(BPDG)))
    if verbose:
        logging.info(f'Xc: {XcPDG}')
    XcID = PROC.add_particle(Particle(FourMomentum(XcE, Xcpx, Xcpy, Xcpz), int(XcPDG)))
    if abs(XcPDG) not in mc_pdg_codes['D']:
        if verbose:
            logging.info('Exciting! An excited D!')
            #logging.info("Xc Dght:", XcD1PDG)
        # add daughters of D* (always 2) or D** (always 3):
        # if this is a narrow D** that decays into a D*, I need to add the two D* daughters as well:
        # this is not necessary for the D daughters as the D has no spin(?). At least D daughters don't
        # change anything.
        #if abs(XcPDG) in [*mc_pdg_codes['D_1'], *mc_pdg_codes['D_2*']] and abs(XcD1PDG) in mc_pdg_codes['D*']:
        #    if verbose:
        #        logging.info("Super exciting! A super excited D decaying into an exciting one!")
    if verbose:
        logging.info(f'Lepton: {lepPDG}')
    lepID = PROC.add_particle(Particle(FourMomentum(lepE, leppx, leppy, leppz), int(lepPDG)))
    # add the remaining B daughters (BD3 is nulep, the others would be photons):
    if verbose:
        logging.info('Add Neutrino')
    BD3ID = PROC.add_particle(Particle(FourMomentum(BD3E, BD3px, BD3py, BD3pz), int(BD3PDG)))

    if abs(XcPDG) not in mc_pdg_codes['D']:
        if verbose:
            logging.info('Exciting! An excited D!')
            logging.info(f'Xc Dght:{XcD1PDG}')
        XcD1ID = PROC.add_particle(Particle(FourMomentum(XcD1E, XcD1px, XcD1py, XcD1pz), int(XcD1PDG)))
        if abs(XcPDG) in [*mc_pdg_codes['D_1'], *mc_pdg_codes['D_2*']] and abs(XcD1PDG) in mc_pdg_codes['D*']:
            if verbose:
                logging.info('Super exciting! A super excited D decaying into an exciting one!')
                logging.info(f'Xc GDght: {XcGD1PDG}')
                logging.info(f'Xc GDght: {XcGD2PDG}')
            XcGD1ID = PROC.add_particle(Particle(FourMomentum(XcGD1E, XcGD1px, XcGD1py, XcGD1pz), int(XcGD1PDG)))
            XcGD2ID = PROC.add_particle(Particle(FourMomentum(XcGD2E, XcGD2px, XcGD2py, XcGD2pz), int(XcGD2PDG)))
            PROC.add_vertex(XcD1ID, list({XcGD1ID, XcGD2ID}))
        if verbose:
            logging.info(f'Xc Dght:{XcD2PDG}')
        XcD2ID = PROC.add_particle(Particle(FourMomentum(XcD2E, XcD2px, XcD2py, XcD2pz), int(XcD2PDG)))
        if not np.isnan(XcD3E):
            if verbose:
                logging.info(f'Xc Dght: {XcD3PDG}')
            XcD3ID = PROC.add_particle(Particle(FourMomentum(XcD3E, XcD3px, XcD3py, XcD3pz), int(XcD3PDG)))
            PROC.add_vertex(XcID, list({XcD1ID, XcD2ID, XcD3ID}))
        else:
            PROC.add_vertex(XcID, list({XcD1ID, XcD2ID}))

    PROC.add_vertex(BID, list({XcID, lepID, BD3ID}))
    # note that sometimes we do have two B mesons, it's up to me if I want to distinguish that
    # (in my case, I took care of that?! But strictly, the tag-side would change as well?!)
    if verbose:
        logging.info('Ready. Setup HAMMER')
    ids = set()
    id_proc = HAMMER.add_process(PROC)  # I think this is only filled if this process is part of hammer
    if (id_proc):
        ids.add(id_proc)
    if (ids):
        HAMMER.process_event()
        # to only account for the shape difference, I'd want to include the overall rate change as well
        # as it is done e.g. here https://stash.desy.de/projects/B2A/repos/sduell_phaseii/browse/scripts/pyHammerModule.py#324-335
        # however, I'd want to think about this again and decide what I'd really want to do
        diffweights = np.array([HAMMER.get_weight(schemename) for schemename in schemename_strings])
        decay_string = pdg_to_string_dict[BPDG] + pdg_to_string_dict[XcPDG] + pdg_to_string_dict[
            lepPDG] + pdg_to_string_dict[BD3PDG]
        if decay_string not in ratedict['denominator'].keys():
            ratedict['denominator'][decay_string] = HAMMER.get_denominator_rate(decay_string)
            logging.info(f'Added denominator rate for {decay_string}: {ratedict["denominator"][decay_string]}')
        decFF_total_rate = ratedict['denominator'][decay_string]
        if np.sum(np.isnan(diffweights)) == 0:
            if verbose:
                logging.info('It works!')
                logging.info(diffweights)
            for i, schemename in enumerate(schemename_strings):
                if decay_string not in ratedict[schemename].keys():
                    ratedict[schemename][decay_string] = HAMMER.get_rate(decay_string, schemename)
                    logging.info(f'Added {schemename} rate for {decay_string}: {ratedict[schemename][decay_string]}')
                newFF_total_rate = ratedict[schemename][decay_string]
                if decFF_total_rate == 0. or newFF_total_rate == 0.:
                    logging.info(
                        f'WARNING! Something in the total rates went wrong?: dec total: {decFF_total_rate}, new total: {newFF_total_rate}, decay_string: {decay_string}'
                    )
                weightdict[schemename].append(decFF_total_rate / newFF_total_rate * diffweights[i])
                #weightdict[schemename + '_decTOT'].append(decFF_total_rate)
                #weightdict[schemename + '_newTOT'].append(newFF_total_rate)
        else:
            if ignore_photons:
                logging.info(
                    f'{np.sum(np.isnan(diffweights))} out of {len(diffweights)} weights are Nan although photons are not considered! The Nans are added.'
                )
                logging.info(f'Failed schemes: {schemename_strings[np.isnan(diffweights)]}')
                for i, schemename in enumerate(schemename_strings):
                    if decay_string not in ratedict[schemename].keys():
                        ratedict[schemename][decay_string] = HAMMER.get_rate(decay_string, schemename)
                        logging.info(
                            f'Added {schemename} rate for {decay_string}: {ratedict[schemename][decay_string]}')
                    newFF_total_rate = ratedict[schemename][decay_string]
                    if decFF_total_rate == 0. or newFF_total_rate == 0.:
                        logging.info(
                            f'WARNING! Something in the total rates went wrong?: dec total: {decFF_total_rate}, new total: {newFF_total_rate}, decay_string: {decay_string}'
                        )
                    weightdict[schemename].append(decFF_total_rate / newFF_total_rate * diffweights[i])
                    #weightdict[schemename + '_decTOT'].append(decFF_total_rate)
                    #weightdict[schemename + '_newTOT'].append(newFF_total_rate)
            else:
                logging.info(
                    f'{np.sum(np.isnan(diffweights))} out of {len(diffweights)} weights are Nan! They are recalculated without photons!'
                )
                logging.info(f'Failed schemes: {schemename_strings[np.isnan(diffweights)]}')
                proc_BtoXcLepNu_event(BE,
                                      Bpx,
                                      Bpy,
                                      Bpz,
                                      BPDG,
                                      XcE,
                                      Xcpx,
                                      Xcpy,
                                      Xcpz,
                                      XcPDG,
                                      XcD1E,
                                      XcD1px,
                                      XcD1py,
                                      XcD1pz,
                                      XcD1PDG,
                                      XcD2E,
                                      XcD2px,
                                      XcD2py,
                                      XcD2pz,
                                      XcD2PDG,
                                      XcD3E,
                                      XcD3px,
                                      XcD3py,
                                      XcD3pz,
                                      XcD3PDG,
                                      XcGD1E,
                                      XcGD1px,
                                      XcGD1py,
                                      XcGD1pz,
                                      XcGD1PDG,
                                      XcGD2E,
                                      XcGD2px,
                                      XcGD2py,
                                      XcGD2pz,
                                      XcGD2PDG,
                                      lepE,
                                      leppx,
                                      leppy,
                                      leppz,
                                      lepPDG,
                                      BD3E,
                                      BD3px,
                                      BD3py,
                                      BD3pz,
                                      BD3PDG,
                                      HAMMER,
                                      weightdict,
                                      ratedict,
                                      ignore_photons=True,
                                      verbose=verbose)
    else:
        for schemename in schemename_strings:
            logging.info(60 * '*', '\nWARNING! Something went wrong as the IDs are not filled?!')
            weightdict[schemename].append(1.)
            #weightdict[schemename + '_decTOT'].append(1.)
            #weightdict[schemename + '_newTOT'].append(1.)
        if verbose:
            logging.info('The process seems not to be part of HAMMER. 1. is used instead.')
    return ratedict


def run_BtoXcLepNu_hammer(df: pd.DataFrame,
                          dict: Dict[str, Union[int, float, np.array]],
                          mode='XcEllNu',
                          verbose=False):
    """
    Functions that initializes HAMMER depending on the mode and runs it.
    :param df: pd.DataFrame that HAMMER is run over.
    :param mode: Mode to specify which models should be added. See initialize_BtoXcLepNu_hammer
                 for more details.
    """
    logging.info(df)
    # df.to_parquet('/afs/desy.de/user/t/tommy/dust/samples/trash/df_pipeline_fail.pq', engine='pyarrow')
    HAMMER = initialize_BtoXcLepNu_hammer(dict, mode=mode, verbose=verbose)
    logging.info('done')
    ff_weight_dict = {}
    tot_rate_dict = {}
    tot_rate_dict['denominator'] = {}
    for schemename in HAMMER.get_ff_scheme_names():
        ff_weight_dict[schemename] = []
        #ff_weight_dict[schemename + '_decTOT'] = []
        #ff_weight_dict[schemename + '_newTOT'] = []
        tot_rate_dict[schemename] = {}

    # if np.any(mask):
    ratedict = df.apply(lambda row: proc_BtoXcLepNu_event(row.Bsig_gen_PDG,
                                                         row.genXBFrameE,
                                                         row.genXBFramePx,
                                                         row.genXBFramePy,
                                                         row.genXBFramePz,
                                                         row.X_gen_PDG,
                                                         row.genXD1BFrameE,
                                                         row.genXD1BFramePx,
                                                         row.genXD1BFramePy,
                                                         row.genXD1BFramePz,
                                                         row.Xdecay0_gen_PDG,
                                                         row.genXD2BFrameE,
                                                         row.genXD2BFramePx,
                                                         row.genXD2BFramePy,
                                                         row.genXD2BFramePz,
                                                         row.Xdecay1_gen_PDG,
                                                         row.genXD3BFrameE,
                                                         row.genXD3BFramePx,
                                                         row.genXD3BFramePy,
                                                         row.genXD3BFramePz,
                                                         row.Xdecay2_gen_PDG,
                                                         row.genXGD1BFrameE,
                                                         row.genXGD1BFramePx,
                                                         row.genXGD1BFramePy,
                                                         row.genXGD1BFramePz,
                                                         row.XGD0_gen_PDG,
                                                         row.genXGD2BFrameE,
                                                         row.genXGD2BFramePx,
                                                         row.genXGD2BFramePy,
                                                         row.genXGD2BFramePz,
                                                         row.XGD1_gen_PDG,
                                                         row.genLepBFrameE,
                                                         row.genLepBFramePx,
                                                         row.genLepBFramePy,
                                                         row.genLepBFramePz,
                                                         row.lep_gen_PDG,
                                                         row.genNeuBFrameE,
                                                         row.genNeuBFramePx,
                                                         row.genNeuBFramePy,
                                                         row.genNeuBFramePz,
                                                         row.neutrino_gen_PDG,
                                                         HAMMER,
                                                         weightdict=ff_weight_dict,
                                                         ratedict=tot_rate_dict,
                                                         ignore_photons=True,
                                                         verbose=verbose),
                       axis=1)
    logging.info(f'Successfully ran over the dataset in the mode {mode}. Now, the weights are added to the dataframe.')
    for schemename in HAMMER.get_ff_scheme_names():
        nans = np.sum(np.isnan(np.array(ff_weight_dict[schemename])))
        if nans > 0:
            if mode in ['XcEllNu', 'XcTauNu']:
                logging.info(
                    f'ERROR in run_BtoXcLepNu_hammer! NaNs were detected but the total size of the concerned decay processes cannot be deduced due to the inclusiveness of the mode {mode}.'
                )
                #raise ValueError(
                #    f'ERROR in run_BtoXcLepNu_hammer! NaNs were detected but the total size of the concerned decay processes cannot be deduced due to the inclusiveness of the mode {mode}.'
                #)
            total_events = len(df)
            logging.info(
                f'In total, {nans} weights with NaN value were detected. They are set to 1e-7 and all other weights are upscaled by {nans}/{total_events}.'
            )
            # NaN + number gives NaN again:
            ff_weight_dict[schemename] = np.array(ff_weight_dict[schemename]) + nans / total_events
            # set NaNs to very small number (zero crashes if you divide by it):
            ff_weight_dict[schemename] = np.nan_to_num(ff_weight_dict[schemename], nan=1e-7)
        df.loc[:, f'FF_weight_{schemename}'] = 1.0
        df.loc[:, f'FF_weight_{schemename}'] = np.array(ff_weight_dict[schemename])
        #df.loc[:, f"__FF_{schemename}_decTOT__"] = np.array(ff_weight_dict[schemename + '_decTOT'])
        #df.loc[:, f"__FF_{schemename}_newTOT__"] = np.array(ff_weight_dict[schemename + '_newTOT'])
    return ratedict


