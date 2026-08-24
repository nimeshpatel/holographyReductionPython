#!/usr/bin/env python3
"""
Nimesh Patel
October 2025
Python version of previous fortran/C code from Holis package.

PREPROCESS.PY - Data preprocessing for holography analysis
Applies calibrations, corrections, and reference plane adjustments.
"""
from pathlib import Path
import math, re
import numpy as np
from scipy import constants

PI = constants.pi
CLIGHT = constants.c  # Speed of light in m/s (exact value: 299792458.0)

PARAM_ORDER_HINTS = [
    ("in_data", "file"),
    ("outamp", "file"),
    ("outphase", "file"),
    ("bore", "file"),
    ("ninp", "int"),
    ("nout", "int"),
    ("frequency", "float"),
    ("samp_itvl", "float"),
    ("dist1", "float"),
    ("dist2", "float"),
    ("ido_taper", "int"),
    ("ido_bore", "int"),
    ("ido_inter", "int"),
    ("adc_of_ph", "float"),
    ("adc_of_am", "float"),
    ("vvm_g_ph", "float"),
    ("dc_g_ph", "float"),
    ("adc_g_ph", "float"),
    ("vvm_g_am", "float"),
    ("dc_g_am", "float"),
    ("adc_g_am", "float"),
]

def parse_preprocess_prm(path: Path):
    text = path.read_text(errors="ignore").splitlines()
    params = {}
    def put(name, value): params[name]=value
    def coerce(val, kind):
        v = val.strip()
        if kind == "file": return v
        if kind == "int": return int(v)
        if kind == "float": return float(v)
        return v
    for line in text:
        if not line.strip() or line.lstrip().startswith('!') or line.strip().startswith('...'):
            continue
        label = line[:40].strip().lower()
        if 'dc' in label and 'ph' in label and 'amplifier' in label:
            print(f"DEBUG: Found DC phase line: '{label}'")

        value_field = line[49:].strip() if len(line) >= 50 else ""
        if not value_field:
            m = re.search(r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*$", line)
            if m: value_field = m.group(1)
            else:
                m2 = re.search(r"([A-Za-z0-9_.+-]+)\s*$", line)
                value_field = m2.group(1) if m2 else ""
        key = None
        lbl = label.replace(' ','').replace('-','').replace('_','')
        if lbl.startswith('input far field amplitude file name') or 'input' in lbl and ('amp' in lbl or 'phase' in lbl or 'rgin' in lbl):
            key='in_data'
        elif lbl.startswith('output far field amplitude file name') or ('output' in lbl and 'amp' in lbl) or lbl.startswith('outamp'):
            key='outamp'
        elif lbl.startswith('output far field phase file name') or ('output' in lbl and 'phase' in lbl) or lbl.startswith('outphase'):
            key='outphase'
        elif lbl.startswith('bore-sight data file name') or 'bore' in lbl:
            key='bore'
        elif lbl.startswith('size ninp of the input data file') or 'ninp' in lbl or ('input' in lbl and 'size' in lbl):
            key='ninp'
        elif lbl.startswith('size nout of the output data file') or 'nout' in lbl or ('output' in lbl and 'size' in lbl):
            key='nout'
        elif lbl.startswith('observing frequency (ghz)') or 'freq' in lbl:
            key='frequency'
        elif lbl.startswith('sampling interval (arcsecs)') or 'samp' in lbl or 'sample' in lbl:
            key='samp_itvl'
        elif 'dist1' in lbl or 'referenceplaneofthemaindish' in lbl or lbl.startswith('reference plane of the main dish'):
            key='dist1'
        elif 'dist2' in lbl or 'referenceplaneofthereferencehorn' in lbl or lbl.startswith('reference plane of the reference horn'):
            key='dist2'
        elif lbl.startswith('do far-field taper') or ('ido' in lbl and 'taper' in lbl):
            key='ido_taper'
        elif lbl.startswith('do bore-sight drift correction') or ('ido' in lbl and 'bore' in lbl):
            key='ido_bore'
        elif lbl.startswith('do interpolation') or ('ido' in lbl and ('inter' in lbl or 'interp' in lbl)):
            key='ido_inter'
        elif 'adcof' in lbl and 'ph' in lbl:
            key='adc_of_ph'
        elif 'adcof' in lbl and ('am' in lbl or 'amp' in lbl):
            key='adc_of_am'

        elif 'vvm' in lbl and 'ph' in lbl:
            key='vvm_g_ph'
        elif 'adc' in lbl and 'gain' in lbl and 'ph' in lbl:
            key='adc_g_ph'
        elif 'dc' in lbl and 'gain' in lbl and 'ph' in lbl:
            key='dc_g_ph'
        elif 'vvm' in lbl and ('am' in lbl or 'amp' in lbl):
            key='vvm_g_am'
        elif 'adc' in lbl and 'gain' in lbl and ('am' in lbl or 'amp' in lbl):
            key='adc_g_am'
        elif 'dc' in lbl and 'gain' in lbl and ('am' in lbl or 'amp' in lbl):
            key='dc_g_am'

        elif lbl.startswith('input units (engineering|counts)'):
            key='units_mode'
        if key:
            kind = next((t for (k,t) in PARAM_ORDER_HINTS if k==key), 'float')
            try: put(key, coerce(value_field, kind))
            except: pass
    # Units override
    if params.get('units_mode','').lower().startswith('engineer'):
        params['units_override'] = 'engineering'
    elif params.get('units_mode','').lower().startswith('count'):
        params['units_override'] = 'counts'

    if all(k in params for k in ('vvm_g_ph','dc_g_ph','adc_g_ph')):
        denom = (params['vvm_g_ph']*params['dc_g_ph']*params['adc_g_ph']) or 1.0
        params['scale_ph'] = 1.0/denom
    if all(k in params for k in ('vvm_g_am','dc_g_am','adc_g_am')):
        denom = (params['vvm_g_am']*params['dc_g_am']*params['adc_g_am']) or 1.0
        params['scale_am'] = 1.0/denom
    params.setdefault('ido_taper',0)
    params.setdefault('ido_bore',0)
    params.setdefault('ido_inter',0)
    params.setdefault('adc_of_ph',0.0)
    params.setdefault('adc_of_am',0.0)
    params.setdefault('dist1',0.0)
    params.setdefault('dist2',0.0)
    required = ['in_data','outamp','outphase','ninp','nout','frequency','samp_itvl']
    missing = [k for k in required if k not in params]
    if missing:
        raise ValueError(f'Missing required parameters in {path.name}: {missing}')
    print(f"DEBUG params keys: {params.keys()}")
    print(f"vvm_g_ph = {params.get('vvm_g_ph')}")
    print(f"dc_g_ph = {params.get('dc_g_ph')}")
    print(f"adc_g_ph = {params.get('adc_g_ph')}")
    params.setdefault('scale_ph', 1.0)
    params.setdefault('scale_am', 1.0)
    return params

def read_raster_and_bore(raster_path: Path, bore_path: Path, ninp: int, use_bore: bool):
    am = [[0.0 for _ in range(ninp)] for _ in range(ninp)]
    ph = [[0.0 for _ in range(ninp)] for _ in range(ninp)]
    with raster_path.open('r') as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith('!') and not ln.strip().startswith('...')]
    idx = 0
    for i in range(ninp):
        for j in range(ninp):
            parts = lines[idx].split()
            if len(parts) >= 4:
                a = float(parts[-2]); p = float(parts[-1])  # last two tokens
            elif len(parts) == 3:
                a = float(parts[1]);  p = float(parts[2])
            else:
                raise ValueError("Bad input line in rgin.dat: " + lines[idx])
            am[i][j] = a
            ph[i][j] = p
            idx += 1

    ambore = [0.0]*(ninp+1)
    phbore = [0.0]*(ninp+1)
    if use_bore and bore_path and Path(bore_path).exists():
        with Path(bore_path).open('r') as f:
            lines_b = [ln.strip() for ln in f if ln.strip() and not ln.strip().startswith('!')]
        for ii in range(ninp+1):
            parts = lines_b[ii].split()
            ambore[ii] = float(parts[-2]); phbore[ii] = float(parts[-1])
    return am, ph, ambore, phbore


def interpol_inplace(am, ph):
    n = len(am); m = len(am[0])
    for i in range(1,n-1):
        for j in range(m):
            if am[i][j] > 38000.0 or am[i][j] < 32500.0:
                am[i][j] = 0.5*(am[i+1][j] + am[i-1][j])
    for i in range(1,n-1):
        for j in range(m):
            if ph[i][j] > 38000.0 or ph[i][j] < 26500.0:
                ph[i][j] = 0.5*(ph[i+1][j] + ph[i-1][j])

def count_conv_inplace(am, ph, ambore, phbore, scale_am, adc_of_am, scale_ph, adc_of_ph, use_bore: bool):
    n = len(am); m = len(am[0])
    for i in range(n):
        for j in range(m):
            am[i][j] = (am[i][j] - adc_of_am)*scale_am
            ph[i][j] = (ph[i][j] - adc_of_ph)*scale_ph
    if use_bore:
        for i in range(len(ambore)):
            ambore[i] = (ambore[i] - adc_of_am)*scale_am
            phbore[i] = (phbore[i] - adc_of_ph)

def apply_bore_inplace(am, ph, ambore, phbore):
    ninp = len(am); ncol = len(am[0])
    for i in range(ninp):
        dr_ph = (phbore[i+1]-phbore[i])/(ninp-1)
        dr_am = (ambore[i+1]-ambore[i])/(ninp-1)
        for j in range(ncol):
            ph[i][j] = ph[i][j] - (dr_ph*j + phbore[i])
            am[i][j] = am[i][j] / (dr_am*j + ambore[i])

def reference_inplace(ph, ninp, dist1, dist2, samp_itvl_arcsec, freq_ghz):
    # Use module-level constants from scipy.constants
    freq = freq_ghz * 1.0e9
    alambda = CLIGHT / freq
    ncent = ninp // 2 + 1
    samp = (samp_itvl_arcsec / 3600.0) * (PI / 180.0)

    delta_d = (dist1 - dist2)
    k = (2.0 * PI / alambda) * delta_d

    # Match Fortran: i outer (row), j inner (column)
    for i in range(ninp):
        for j in range(ninp):
            # Fortran 1-based: x = j-ncent, y = i-ncent where j,i start at 1
            x = float((j + 1) - ncent)
            y = float((i + 1) - ncent)
            radius = math.sqrt(x * x + y * y)
            theta = radius * samp
            ph[i][j] += k * (1.0 - math.cos(theta))

def sizeconv(am, ph, ninp, nout):
    am1 = [[0.0 for _ in range(nout)] for _ in range(nout)]
    ph1 = [[0.0 for _ in range(nout)] for _ in range(nout)]
    ndiff = (nout//2) - (ninp//2)
    n1 = ndiff
    n2 = nout - ndiff
    for i in range(n1, n2):
        for j in range(n1, n2):
            am1[i][j] = am[i-n1][j-n1]
            ph1[i][j] = ph[i-n1][j-n1]
    return am1, ph1


def autodetect_units(amps, phases):
    # If amplitudes are small (<1e2) and phases lie within [-pi, pi], assume engineering units (V, rad).
    import math
    flat_amp = [v for row in amps for v in row]
    flat_ph = [v for row in phases for v in row]
    if not flat_amp or not flat_ph:
        return "counts"
    amax = max(abs(v) for v in flat_amp)
    pmax = max(abs(v) for v in flat_ph)
    if amax < 100.0 and pmax <= math.pi*1.05:
        return "engineering"
    return "counts"

def write_data(am, ph, outamp_path: Path, outphase_path: Path):
    """Write data in the same order it was computed (row-major)."""
    with Path(outamp_path).open('w') as fa, Path(outphase_path).open('w') as fp:
        n = len(am)
        m = len(am[0])
        for i in range(n):          # outer = elevation (same as compute)
            for j in range(m):      # inner = azimuth
                fa.write(f"{am[i][j]:.16e}\n")
                fp.write(f"{ph[i][j]:.16e}\n")


def preprocess_main(prm_path: Path):
    prm = parse_preprocess_prm(prm_path)
    print(f"scale_ph = {prm.get('scale_ph')}, scale_am = {prm.get('scale_am')}")  # DEBUG
    print(f"adc_of_ph = {prm.get('adc_of_ph')}, adc_of_am = {prm.get('adc_of_am')}")  # DEBUG
    ninp = int(prm['ninp']); nout = int(prm['nout'])
    raster_path = Path(prm['in_data'])
    bore_path = Path(prm.get('bore',''))
    outamp_path = Path(prm['outamp'])
    outphase_path = Path(prm['outphase'])
    use_bore = int(prm.get('ido_bore',0)) == 1

    am, ph, ambore, phbore = read_raster_and_bore(raster_path, bore_path, ninp, use_bore)

    print(f"After reading: ph[0][0] = {ph[0][0]}")  # DEBUG

    # Auto-detect units
    units = prm.get('units_override') or autodetect_units(am, ph)
    print(f"Units detected: {units}")  # DEBUG

    if units == 'counts' and int(prm.get('ido_inter',0)) == 1:
        interpol_inplace(am, ph)
    count_conv_inplace(am, ph, ambore, phbore,
                           prm['scale_am'], prm['adc_of_am'],
                           prm['scale_ph'], prm['adc_of_ph'],
                           use_bore)
    
    print(f"After count_conv: ph[0][0] = {ph[0][0]}")  # DEBUG

    if use_bore:
        apply_bore_inplace(am, ph, ambore, phbore)

    print(f"After bore: ph[0][0] = {ph[0][0]}")  # DEBUG

    reference_inplace(ph, ninp, prm.get('dist1',0.0), prm.get('dist2',0.0),
                      prm['samp_itvl'], prm['frequency'])

    print(f"After reference: ph[0][0] = {ph[0][0]}")  # DEBUG

    am1, ph1 = sizeconv(am, ph, ninp, nout)
    print(f"After sizeconv: ph1[0][0] = {ph1[0][0]}")  # DEBUG
    write_data(am1, ph1, outamp_path, outphase_path)

if __name__ == '__main__':
    preprocess_main(Path('preprocess.prm'))
