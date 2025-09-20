import pandas as pd
import math
import numpy as np

def rank_n(variable, peso, n_threshold=3, best_threshold=0.3, middle_threshold=0.7):
    """
    Asigna un puntaje escalado en función de umbrales.
    Un valor menor es considerado más favorable.

    Args:
        variable (float): valor de la propiedad.
        peso (float): peso asignado al criterio.
        n_threshold (int | float): número de niveles de puntuación.
        best_threshold (float): umbral óptimo.
        middle_threshold (float): umbral intermedio.

    Returns:
        float: puntaje parcial ajustado al peso.
    """
    if variable is None:
        return 0.0

    if n_threshold == 3:
        if variable <= best_threshold:
            return peso * 1.0
        elif variable <= middle_threshold:
            return peso * 0.5
        else:
            return peso * 0.3

    elif n_threshold == 2:
        if variable <= best_threshold:
            return peso * 1.0
        elif variable <= middle_threshold:
            return peso * 0.5

    elif n_threshold == 1.5:
        return peso * 1.0 if variable <= best_threshold else peso * 0.5

    elif n_threshold == 1:
        return peso * 1.0 if variable <= best_threshold else 0.0

    return 0.0

def calcular_puntaje_admet(row, pesos=None, return_breakdown=False,  verbose=False, FQ=True, MedChem=True, Abs=True, Dis=True, Met=True, Exc=True, Tox=True,):
    """
    Calcula un puntaje compuesto de ADMET para priorizar compuestos.

    Args:
        row (pd.Series): Fila del DataFrame con propiedades ADMET.
        pesos (dict, optional): Pesos personalizados para cada criterio.
        verbose (bool): Si True, imprime los detalles del puntaje.
        return_breakdown (bool): Si True, devuelve desglose de contribuciones.

    Returns:
        float: Puntaje total normalizado (0-1).
    """
    DEFAULT_WEIGHTS = {
#        #FQ
        'stereo': 1.0, #FQ
        'logD': 1.0, #ADMET_LAB3
        'logP':1.0,
        'free_energy_hyd': 1.0,

#        #Medicinal chemistry
        'qed': 0.5, 'SAscore':1.0,
        'ghose': 1.0, 'lead': 1.0,'brenk': 1.0, 'bioavailability': 0.5, #SwissAdme
        'alarm_nmr': 1.0,'coloidal_ag': 1.0, #LAB3
        'reactive': 0.5, 'green_fluor':0.5, 'other_interference':0.5, #Lab3

#        #Absorción
        'caco2': 1.5, 'pampa':1.5,
        'mdck': 1.0, #lab3
        'pgp_subs': 1.5, 'hia': 0.5,#pkcsm

#       #Distribución
        'logvds': 1.0,'fu': 0.5,'cns_permeability': 0.5, #pkcsm
        'mrp1': 1.0, 'bsep':1.0, 'oatp1b3':1.0, 'oatp1b1':1.0, #Lab3

#        #Metabolismo
        'CYP2C9_inh':0.5, 'CYP2C8_inh': 1.0,#lab3
        'CYP2C19_inh':1.0, 'CYP2C19_sub': 1.0,#lab3
        'CYP2D6_inh': 1.5, 'CYP2D6_sub': 1.5,
        'CYP3A4_inh': 1.5,'CYP3A4_sub': 1.5,
        'CYP2B6_inh': 1.0, 'hlm_stab': 1.0, #Lab3

#        #Excreción; no se considera t05
        'oct2_subs':1.0, 'cl_total':0.5, #pkcsm
        'cl_plasma':0.5, #lab3

#        #Toxicidad
        'neuro_tox':1.0, 'oto_tox':1.0,'hema_tox': 1.0, 'gen_tox': 1.0,'carcino_tox': 1.0,  #Lab3
        'ames': 1.5, 'herg': 0.5,'herg_10um': 0.5,'hepa_tox': 1.0, 'dili': 1.0,'LC50FM':1.0,
        'A549':1.0, 'HEK293':0.5, 'ROA':1.0, 'FDAMDD':1.0, 'IGC50':1.0, 'LC50DM':1.0, #Lab3
        'nr_ar':1.0, 'nr_ar_lbd':0.5, 'nr_er':1.0, 'nr_arom':1.0, 'nr_ahr': 0.5, 'sr_are':1.0, #lab3
        'sr_atad5': 0.5, 'sr_mmp':0.5, 'sr_p53':1.0, 'sure_chembl':1.0, #Lab3
    }
#
    pesos = pesos or DEFAULT_WEIGHTS
    if not isinstance(row, pd.Series):
        raise ValueError("Entrada esperada: una fila del DataFrame (pd.Series).")

    score = 0.0
    max_score = 0.0
    breakdown = {}

    def add_score(key, value, use_rank=False, 
                  FQ=True, MedChem=True, Abs=True, Dis=True, Met=True, Exc=True, Tox=True,
                  **kwargs):
        '''
        Función interna que simplifica la suma del puntaje parcial para cada propiedad
        ADMET.
        key: nombre del criterio (propiedad)
        value: valor de la propiedad para esa fila
        use_rank: indica si se debe de aplicar la función rank_n (True) o solo se
          usa el peso completo (False)
        **kwargs: argumentos adicionales para la función rank_n
        '''
        nonlocal score, max_score #Funciones no locales a add_score, pertenen al mundo externo
        peso = pesos.get(key, 0.0) #Peso asignado al criterio
        if peso == 0.0  or value is None:
            return
        if value is not None:
            contrib = rank_n(value, peso, **kwargs) if use_rank else peso #Si es True, calcula contribución
            score += contrib
            max_score += peso
            breakdown[key] = contrib  # << registrar contribución

            if verbose:
                print(f"Propiedad: {key}, Valor: {value}, Peso: {peso}, Contribución: {contrib:.2f}")


    try:


#        ''' ===== FQ ===== # : '''
      if FQ: 
        
        add_score('stereo', row.get('nStereo', row.get('stereo_centers')), use_rank=True, n_threshold=1, best_threshold=2) #Lab3, AI

        logd_value = row.get('logD') #lab3:[ 'stereo': 1.0, ]#logD logP
        peso_logd = pesos.get('logD', 0.0)
        if logd_value is not None and peso_logd > 0:
            if 1 <= logd_value <= 3:
                contrib = peso_logd
            else:
                contrib = peso_logd * 0.5
            score += contrib
            max_score += peso_logd
            breakdown['logD'] = contrib

        logp_value = row.get('WLOGP') or row.get('LOGP') or row.get('logP')
        peso_logp = pesos.get('logP', 0.0)

        if logp_value is not None and peso_logp > 0:
            if 0.95 <= logp_value <= 2.2:
                contrib = peso_logp
            elif logp_value >= 0:
                contrib = peso_logp * 0.5
            else:
                contrib = 0.0
            score += contrib
            max_score += peso_logp
            breakdown['logP'] = contrib

        free_energy_hyd = row.get('HydrationFreeEnergy_FreeSolv')
        peso_hyd = pesos.get('free_energy_hyd', 0.0)

        if free_energy_hyd is not None and peso_hyd > 0:
            if free_energy_hyd >= -7:
                contrib = peso_hyd
            else:
                contrib = peso_hyd * 0.5

            score += contrib
            max_score += peso_hyd
            breakdown['free_energy_hyd'] = contrib

#        ''' ===== MEDICINAL CHEMISTRY ===== # : '''
      if MedChem: 
        qed = row.get('QED') # Quantitative Estimate of Drug-likeness // Lab3, AI
        peso_qed = pesos.get('qed', 0.0)
        if qed is not None and peso_qed > 0:
            if qed >= 0.67:
                contrib = peso_qed
            elif qed > 0.49:
                contrib = peso_qed * 0.5
            else:
                contrib = 0.0
            score += contrib 
            max_score += peso_qed
            breakdown['qed'] = contrib


        add_score('SAscore', row.get('Synthetic Accessibility', row.get('Synth')), use_rank=True, n_threshold=2, best_threshold=4.25, middle_threshold=6)

        add_score('ghose', row.get('Ghose #violations'), use_rank=True, n_threshold=2, best_threshold=0, middle_threshold=1)# Swiss

        add_score('lead', row.get('Leadlikeness #violations'), use_rank=True, n_threshold=2, best_threshold=0, middle_threshold=1)# Swiss

        add_score('brenk', row.get('Brenk #alerts'), use_rank=True, n_threshold=2, best_threshold=0, middle_threshold=1)#Alertas estructurales // # Swiss

        bioavailability = row.get('Bioavailability_Ma', row.get('Bioavailability Score'))
        peso_bio = pesos.get('bioavailability', 0.0)
        if bioavailability is not None and peso_bio > 0:
            contrib = peso_bio if bioavailability >= 0.55 else peso_bio*0.5
            score += contrib
            max_score += peso_bio
            breakdown['bioavailability'] = contrib

        alarm_nmr = row.get('Alarm_NMR') #compuestos reactivos a tiol
        peso_alarm_nmr = pesos.get('alarm_nmr', 0.0)
        if alarm_nmr is not None and peso_alarm_nmr > 0:
              if (isinstance(alarm_nmr, list) and alarm_nmr == ['-']) or str(alarm_nmr) == "['-']":
                  contrib = peso_alarm_nmr
              else:
                  contrib = peso_alarm_nmr * 0.5
              score += contrib
              max_score += peso_alarm_nmr
              breakdown['alarm_nmr'] = contrib

        add_score('coloidal_ag', row.get('Aggregators'), 
                  use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7)# Lab3

        add_score('reactive', row.get('Reactive'), 
                  use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7)

        add_score('green_fluor', row.get('Green_fluorescence'), 
                  use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7)

        add_score('other_interference', row.get('Other_assay_interference'), 
                  use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7)


#        ''' ===== ABSORCIÓN ===== # : '''
      if Abs: 
        caco2_papp =row.get('Caco2 permeability')
        caco2 = row.get('Caco2_Wang')  or row.get('caco2') # Permeabilidad Caco-2 (logPapp)  >= -5.15 //  Lab3 AI
        peso_caco2 = pesos.get('caco2', 0.0)

        if caco2 is not None and peso_caco2 > 0:
          if caco2 >= -5.15: #Mayor valor en refs
              contrib = peso_caco2
          else: 
              contrib =  peso_caco2 * 0.5
          score += contrib
          max_score += peso_caco2
          breakdown['caco2'] =contrib

        if caco2_papp is not None and peso_caco2 > 0:
          if caco2_papp >= 0.9:
            contrib = peso_caco2
          else:
            contrib = peso_caco2 * 0.5
          score += contrib
          max_score += peso_caco2
          breakdown['caco2'] =contrib

        #PAMPA
        pampa = row.get('PAMPA_NCATS') #Lab3
        peso_pampa = pesos.get('pampa', 0.0) 
        if pampa is not None and peso_pampa > 0:
          if  pampa >= 0.7:
            contrib = peso_pampa
          elif pampa >= 0.3:
            contrib = peso_pampa * 0.5
          else:
            contrib = 0.0
          score += contrib
          max_score += peso_pampa
          breakdown['pampa'] = contrib

        add_score('pampa', row.get('PAMPA') , use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7)
        #MDCK
        mdck = row.get('MDCK') # Permeabilidad  >= 2E-6// Lab3
        peso_mdck = pesos.get('mdck', 0.0)
        if mdck is not None and peso_mdck > 0:
          if mdck >= -5.26 and mdck<= 5.17:
            contrib = peso_mdck
          else:
            contrib = peso_mdck * 0.5
          score += contrib
          max_score += peso_mdck
          breakdown['mdck'] = contrib

        add_score('pgp_subs', row.get('pgp_sub'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3

        pgp_substrate= row.get('Pgp substrate',row.get('P-glycoprotein substrate')) #Swiss, PKCSM
        peso_pgp = pesos.get('pgp_subs', 0.0)
        
        if pgp_substrate is not None and peso_pgp > 0:
          if pgp_substrate == 'No':
            contrib = peso_pgp
          else:
            contrib = peso_pgp * 0.5
          score += contrib
          max_score += peso_pgp
          breakdown['pgp_subs'] = contrib

        hia = row.get('Intestinal absorption (human)')
        peso_hia = pesos.get('hia', 0.0)
        if hia is not None and peso_hia > 0:
          if hia >= 30 and hia <= 70:
            contrib = peso_hia
          else:
            contrib = peso_hia * 0.5
          score += contrib
          max_score += peso_hia
          breakdown['hia'] = contrib

#        ''' ===== DISTRIBUCIÓN ===== # : '''
      if Dis: 
        logvds = row.get('VDss (human)') #pkcsm
        peso_logvds = pesos.get('logvds', 0.0)
        if logvds is not None and peso_logvds > 0:
          if logvds >= 0.45 and logvds <= 0.75:
            contrib = peso_logvds
          else:
            contrib = peso_logvds * 0.5
          score += contrib
          max_score += peso_logvds
          breakdown['logvds'] = contrib

        logvds_lombardo = row.get('VDss_Lombardo' ) #AI
        if logvds_lombardo is not None and peso_logvds > 0:
          if  logvds_lombardo >= 2 and logvds_lombardo <= 5:
            contrib = peso_logvds
          else:
            contrib = peso_logvds * 0.5
          score += contrib
          max_score += peso_logvds
          breakdown['logvds'] = contrib

        logvds_lab3 = row.get('logVDss') #Lab3
        if logvds_lab3 is not None and peso_logvds > 0:
          if  logvds_lab3 >= 0.04:
            contrib = peso_logvds
          else:
            contrib = peso_logvds * 0.5
          score += contrib
          max_score += peso_logvds
          breakdown['logvds'] = contrib

        #Fracción no unida (FU):
        fu = row.get('Fraction unbound (human)')
        peso_fu = pesos.get('fu', 0.0)
        if fu is not None and peso_fu > 0:
          if fu >= 0.1 and fu <=0.5:
            contrib = peso_fu
          else:
            contrib = peso_fu * 0.5
          score += contrib
          max_score += peso_fu
          breakdown['fu'] = contrib

        fu_lab3 = row.get('Fu')
        if fu_lab3 is not None and peso_fu > 0:
          if fu_lab3 >= 58.3 and fu_lab3 <= 62.0:
            contrib = peso_fu
          else:
            contrib = peso_fu * 0.5
          score += contrib
          max_score += peso_fu
          breakdown['fu'] = contrib

        '''
        #NO SE CONSIDERA EN EL RANKING DEBIDO A QUE NO TENEMOS EL CONTEXTO DE USO --->¿QUEREMOS QUE ATRAVIESE O NO LA BBB?
        ## NO SE CONSIDERA PPB NI BCRP EN EL RANKING.

        add_score('bbb_perm', row.get('BBB', row.get('BBB_Martins', None)), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)
        bbb_class = row.get('BBB permeant', None)
        if bbb_class is not None:
          if bbb_class == 'No':
            score += pesos['bbb_perm']
            max_score += pesos['bbb_perm']
            breakdown['bbb_perm'] = pesos['bbb_perm']
          else:
            score += pesos['bbb_perm'] * 0.5
            max_score += pesos['bbb_perm']
            breakdown['bbb_perm'] = pesos['bbb_perm'] * 0.5'''

        cns_perm = row.get('CNS permeability') #pkcsm
        peso_cns_perm = pesos.get('cns_permeability', 0.0)
        if cns_perm is not None and peso_cns_perm > 0:
          if cns_perm <= -2:
            contrib = peso_cns_perm
          else:
            contrib = peso_cns_perm * 0.5
          score += contrib
          max_score += peso_cns_perm
          breakdown['cns_permeability'] = contrib

        add_score('mrp1', row.get('MRP1'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('bsep', row.get('BSEP'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('oatp1b3', row.get('OATP1B3'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('oatp1b1', row.get('OATP1B1'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3


#        ''' ===== METABOLISMO ===== # : '''
      if Met: 

        add_score('CYP2C9_inh', row.get('CYP2C9-inh'), use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('CYP2C8_inh', row.get('CYP2C8-inh'), use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('CYP2C19_inh', row.get('CYP2C19-inh'), use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('CYP2C19_sub', row.get('CYP2C19-sub'), use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7) #Lab3
        
        cyp2c19_inh = row.get('CYP2C19 inhibitor')
        peso_cyp2c19_inh = pesos.get('CYP2C19_inh', 0.0)
        if cyp2c19_inh is not None and peso_cyp2c19_inh > 0:
          if  cyp2c19_inh == 'No':
            contrib = peso_cyp2c19_inh
          else:
            contrib = peso_cyp2c19_inh * 0.5
          score += contrib
          max_score += peso_cyp2c19_inh
          breakdown['CYP2C19_inh'] = contrib

        add_score('CYP2D6_inh', row.get('CYP2D6-inh'), use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7) #Lab3
        add_score('CYP2D6_sub', row.get('CYP2D6-sub', row.get('CYP2D6_Substrate_CarbonMangels')), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)

        add_score('CYP3A4_inh', row.get('CYP3A4-inh', row.get('CYP3A4_Veith')), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3

        add_score('CYP3A4_sub', row.get('CYP3A4-sub', row.get('CYP3A4_Substrate_CarbonMangels')), use_rank=True, n_threshold=2, best_threshold=0.3, middle_threshold=0.7)

        cyp3a4_subs = row.get('CYP3A4 substrate')
        peso_cyp3a4_subs = pesos.get('CYP3A4_sub', 0.0)
        if cyp3a4_subs is not None and peso_cyp3a4_subs > 0:
          if cyp3a4_subs == 'No':
            contrib = peso_cyp3a4_subs
          else:
            contrib = peso_cyp3a4_subs * 0.5
          score += contrib
          max_score += peso_cyp3a4_subs
          breakdown['CYP3A4_sub'] = contrib

        add_score('CYP2B6_inh', row.get('CYP2B6-inh'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3

        add_score('hlm_stab', row.get('LM-human'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7) #Lab3

#        ''' ===== EXCRECIÓN ===== # : '''
      if Exc: 
        oct2_subs = row.get('Renal OCT2 substrate' )
        peso_oct2_subs = pesos.get('oct2_subs', 0.0)
        if oct2_subs is not None and peso_oct2_subs > 0:
          if  oct2_subs == 'No':
            contrib = peso_oct2_subs
          else:
            contrib = peso_oct2_subs * 0.5
          score += contrib
          max_score += peso_oct2_subs
          breakdown['oct2_subs'] = contrib

        cl_total = row.get('Total Clearance')
        peso_cl_total = pesos.get('cl_total', 0.0)
        if cl_total is not None and peso_cl_total > 0:
          if cl_total >= 1.2 and cl_total <= 1.3:
            contrib = peso_cl_total
          else:
            contrib = peso_cl_total * 0.5
          score += contrib
          max_score += peso_cl_total
          breakdown['cl_total'] = contrib

        cl_plasma = row.get('cl-plasma')
        peso_cl_plasma = pesos.get('cl_plasma', 0.0)
        if cl_plasma is not None and peso_cl_plasma > 0:
          if cl_plasma >= 7 and cl_plasma <= 8:
            contrib = peso_cl_plasma
          else:
            contrib = peso_cl_plasma * 0.5
          score += contrib
          max_score += peso_cl_plasma
          breakdown['cl_plasma'] = contrib

#        ''' ===== TOXICIDAD ===== # : '''
        #No considero Lab3: FAD-Drugs4 rule, Nefrotoxicity (no considerado por rango (0.98-1.0)), factor de bioconcetracion,
            #pkcsm:  dosis max tol,
            #AI: LD50 ratas
      if Tox: 

        add_score('resp_tox', row.get('Respiratory'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('neuro_tox', row.get('Neurotoxicity-DI'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('oto_tox', row.get('Ototoxicity'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('hema_tox', row.get('Hematotoxicity'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('gen_tox', row.get('Genotoxicity'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3

        add_score('ames', row.get('Ames'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('ames', row.get('AMES'), use_rank=True, n_threshold=3, best_threshold=0.1, middle_threshold=0.2)#IA
        ames_tox = row.get('AMES toxicity')
        peso_ames = pesos.get('ames', 0.0)
        if ames_tox is not None and peso_ames > 0:
          if ames_tox == 'No':
            contrib = peso_ames
          else:
            contrib = peso_ames * 0.5
          score += contrib
          max_score += peso_ames
          breakdown['ames'] = contrib

        add_score('carcino_tox', row.get('Carcinogenicity'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('herg', row.get('hERG'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)
        add_score('herg_10um', row.get('hERG-10um'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('hepa_tox', row.get('H-HT'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('dili', row.get('DILI'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.8)#Lab3
        add_score('dili', row.get('DILI_AI'), use_rank=True, n_threshold=3, best_threshold=0.05, middle_threshold=0.10)#AI ---- CAMBIAR NOMBRE DE DILI A DILI_AI
        dili_pkcsm = row.get('Hepatotoxicity')
        peso_dili_pkcsm = pesos.get('dili', 0.0)
        if dili_pkcsm is not None and peso_dili_pkcsm > 0:
          if dili_pkcsm == 'No':
            contrib = peso_dili_pkcsm
          else:
            contrib = peso_dili_pkcsm * 0.5
          score += contrib
          max_score += peso_dili_pkcsm
          breakdown['dili'] = contrib

        add_score('A549', row.get('A549'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('HEK293', row.get('HEK293'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('ROA', row.get('ROA'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('FDAMDD', row.get('FDAMDD'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        
        igc50_value = row.get('IGC50')
        peso_igc50 = pesos.get('IGC50', 0.0)
        if igc50_value is not None and peso_igc50 > 0:
          if  igc50_value <= 3.5:
            contrib = peso_igc50
          else:
            contrib = peso_igc50 * 0.5
          score += contrib
          max_score += peso_igc50
          breakdown['IGC50'] = contrib

        add_score('IGC50', row.get('<i>T.Pyriformis</i> toxicity'), use_rank=True, n_threshold=1, best_threshold=0.3)#pkcsm
        add_score('LC50FM', row.get('LC50FM'), use_rank=True, n_threshold=3, best_threshold=3, middle_threshold=4.3)#Lab3
        
        lc50fm_value = row.get('Minnow toxicity')
        peso_lc50fm = pesos.get('LC50FM', 0.0)
        if lc50fm_value is not None and peso_lc50fm > 0:
          if  lc50fm_value >=2.06:
            contrib = peso_lc50fm
          else:
            contrib = peso_lc50fm * 0.5
          score += contrib
          max_score += peso_lc50fm
          breakdown['LC50FM'] = contrib

        add_score('LC50DM', row.get('LC50DM'), use_rank=True, n_threshold=3, best_threshold=4.5, middle_threshold=4.7)#Lab3
        add_score('nr_ar', row.get('NR-AR'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('nr_ar_lbd', row.get('NR-AR-LBD'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('nr_er', row.get('NR-ER'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('nr_arom', row.get('NR-Aromatase'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('nr_ahr', row.get('NR-AhR'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('sr_are', row.get('SR-ARE'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('sr_atad5', row.get('SR-ATAD5'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('sr_mmp', row.get('SR-MMP'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3
        add_score('sr_p53', row.get('SR-p53'), use_rank=True, n_threshold=3, best_threshold=0.3, middle_threshold=0.7)#Lab3

        sure_chembl = row.get('SureChEMBL')
        peso_sure_chembl = pesos.get('sure_chembl', 0.0)
        if sure_chembl is not None and peso_sure_chembl > 0:
          if (isinstance(sure_chembl, list) and sure_chembl == ['-']) or str(sure_chembl) == "['-']":

            contrib = peso_sure_chembl
          else:
            contrib = peso_sure_chembl * 0.5
          score += contrib
          max_score += peso_sure_chembl
          breakdown['sure_chembl'] = contrib                   
#
        if verbose:
            smiles = row.get('smiles', 'N/A')
            print(f"Mol: {smiles}, Puntaje: {score:.2f}/{max_score:.2f}")

    except Exception as e:
        if verbose:
            print(f"Error procesando fila: {e}")
        return 0.0

    norm_score = round(score / max_score, 4) if max_score > 0 else 0.0
    if return_breakdown:
        return score, norm_score, breakdown   # 👈 SIEMPRE una tupla
    return norm_score


#Ejemplo de uso:
dfs=[df1, df2, df3]
for df in dfs:
    results = df.apply(calcular_puntaje_admet, axis=1, return_breakdown=True, 
                       FQ=False, MedChem=False, Abs=False, Dis=False, Met=False, Exc=False, Tox=True)
    df['ADMET_Score'] = results.apply(lambda x: x[0])       # el score normalizado
    df['ADMET_Norm_Score'] = results.apply(lambda x: x[1])       # el score normalizado
    df['ADMET_Breakdown'] = results.apply(lambda x: x[2])   # el detalle en dict
    df['ADMET_Rank'] = df['ADMET_Score'].rank(ascending=False, method='min')
