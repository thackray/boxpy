froms = {
    'atm':{
        'atm': 0.,
        'land_arm': 0.05,
        'land_fast': 0.4,
        'land_slow': 0.08,
        'nor_at': 0.7*0.037,
        'deep_nor_at': 0.,
        'sfc_at': 0.7*0.17,
        'int_at': 0.7*0.072,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.7*0.0053,
        'deep_med': 0.,
        'sou_oc': 0.7*0.027,
        'deep_sou_oc': 0.,
        'int_pac': 0.7*0.176,
        'deep_pac': 0.,
        'sfc_pac': 0.7*0.42,
        'nor_pac': 0.7*0.09,
        },
    
    'land_arm':{
        'atm': 2.e-4,
        'land_arm': 0.,
        'land_fast': 8.e-5,
        'land_slow': 0.,
        'nor_at': 3.e-5*0.02,
        'deep_nor_at': 0.,
        'sfc_at': 3.e-5*0.36,
        'int_at': 3.e-5*0.02,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 3.e-5*0.04,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 3.e-5*0.26,
        'nor_pac': 3.e-5*0.30,
        },
    
    'land_fast':{
        'atm': 0.2,
        'land_arm': 1.e-3,
        'land_fast': 0.,
        'land_slow': 0.03,
        'nor_at': 0.04*0.02,
        'deep_nor_at': 0.,
        'sfc_at': 0.04*0.36,
        'int_at': 0.04*0.02,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.04*0.04,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.04*0.26,
        'nor_pac': 0.04*0.30,
        },
    
    'land_slow':{
        'atm': 7.e-3,
        'land_arm': 1.e-5,
        'land_fast': 6.e-3,
        'land_slow': 0.,
        'nor_at': 3.e-4*0.02,
        'deep_nor_at': 0.,
        'sfc_at': 3.e-4*0.36,
        'int_at': 3.e-4*0.02,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 3.e-4*0.04,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 3.e-4*0.26,
        'nor_pac': 3.e-4*0.30,
        },
    
    'nor_at':{
        'atm': 6.e-3,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.02,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_nor_at':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': -3.e-4,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.02,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'sfc_at':{
        'atm': 0.1,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.03,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.1,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 1.e-3,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'int_at':{
        'atm': 8.e-3,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 2.e-3,
        'int_at': 0.,
        'deep_at': 8.e-3,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_at':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 3.e-3,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 5.e-3,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'bot_at':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 4.e-3,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': -7.e-3,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'med_sea':{
        'atm': 0.2,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.05,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_med':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 8.e-3,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': -1.e-3,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'sou_oc':{
        'atm': 1.e-3,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.02,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.02,
        'int_pac': 0.,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_sou_oc':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 6.e-3,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.01,
        'deep_sou_oc': -2.e-4,
        'int_pac': 0.,
        'deep_pac': 0.02,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'int_pac':{
        'atm': 0.01,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 7.e-4,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 4.e-3,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'deep_pac':{
        'atm': 0.,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 3.e-4,
        'deep_pac': -9.e-4,
        'sfc_pac': 0.,
        'nor_pac': 7.e-4,
        },
    
    'sfc_pac':{
        'atm': 0.1,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 0.,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 9.e-3,
        'deep_sou_oc': 0.,
        'int_pac': 0.07,
        'deep_pac': 0.,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    
    'nor_pac':{
        'atm': 0.02,
        'land_arm': 0.,
        'land_fast': 0.,
        'land_slow': 0.,
        'nor_at': 0.,
        'deep_nor_at': 0.,
        'sfc_at': 9.e-3,
        'int_at': 0.,
        'deep_at': 0.,
        'bot_at': 0.,
        'med_sea': 0.,
        'deep_med': 0.,
        'sou_oc': 0.,
        'deep_sou_oc': 0.,
        'int_pac': 0.,
        'deep_pac': 4.e-3,
        'sfc_pac': 0.,
        'nor_pac': 0.,
        },
    }

