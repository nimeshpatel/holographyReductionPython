# gltDishMap.py Revisions for Holography Map Consistency

## Changes Made (2025-10-28)

### Issue Identified
The Python-generated holography maps (`holomapPython.pdf`) showed differences compared to FORTRAN PGPLOT-generated maps:
1. **Panel boundaries misaligned** - surface errors not aligned with panel edges
2. **Masking differences** - quadrupod legs appeared thicker and central hole larger in Python version

### Root Causes Found

1. **Incorrect pixel spacing**:
   - Old code used hardcoded `140.0 mm`
   - Correct value: `124.625 mm` (calculated from holography parameters)
   - Calculation: `(dprim/rate) * 1000 / grid_size = (12.0/0.7521) * 1000 / 128`

2. **Simple circular mask instead of detailed mask**:
   - Old code: Simple `r > 6000 mm` circular mask
   - New code: Reads actual mask from `mask32.dat` including:
     - Outer edge (11.9m diameter)
     - Subreflector hole (0.75m diameter)
     - Quadrupod legs (0.0375m half-width)

3. **Hardcoded parameters**:
   - Old code hardcoded rate, dprim, frequency
   - New code reads from `withphase_aber.prm`

### Changes Implemented

#### 1. Added Parameter File Reader
```python
def read_holography_params(prm_file='withphase_aber.prm'):
```
- Reads Nyquist rate, primary diameter, frequency, grid size
- Falls back to defaults if file not found
- Ensures consistency with holography processing

#### 2. Corrected Pixel Spacing Calculation
```python
fov_meters = dprim / rate
pixel_spacing_mm = (fov_meters * 1000.0) / grid_size
```
- Now calculates from actual parameters: **124.625 mm** instead of 140 mm
- This is a **11% scale correction** - explains panel misalignment!

#### 3. Use Actual Holography Mask
```python
mask_data = numpy.loadtxt(mask_file, dtype=int)
mask_grid = mask_data.reshape((grid_size, grid_size))
mask = (mask_grid == 0)  # Convert to matplotlib mask convention
```
- Reads `mask32.dat` from holography processing
- Shows exact quadrupod legs (vertical/horizontal cross pattern)
- Shows correct subreflector central hole size
- Falls back to simple circular mask if file not available

#### 4. Corrected Error Scaling
```python
wavelength_microns = 299792458.0 / (freq_ghz * 1e9) * 1e6
error_scaling = wavelength_microns / (4.0 * numpy.pi)
```
- Computes from actual frequency (94.5 GHz → 3171 μm wavelength)
- Surface error = λ/(4π) × phase error [radians]
- Result: **252.45 μm/radian** (was 504.90)

#### 5. New Command-Line Arguments
```bash
--mask-file MASK     Mask file from holography (default: mask32.dat)
--prm-file PRM       Parameter file (default: withphase_aber.prm)
```

### Expected Improvements

After these changes, Python-generated maps should show:

1. ✅ **Panel boundaries aligned** with surface errors due to correct pixel spacing
2. ✅ **Thinner quadrupod legs** (37.5mm half-width) matching FORTRAN
3. ✅ **Smaller central hole** (375mm radius) matching FORTRAN
4. ✅ **Correct scale** - 15.95m field of view, not 17.92m
5. ✅ **Consistent with holography processing** - uses same mask and parameters

### Usage

**Basic usage (uses default files):**
```bash
python gltDishMap.py Epr.dat --output holography_map.pdf --vmin -300 --vmax 300 \
       --cmap seismic --label "Surface Error [μm]"
```

**With custom files:**
```bash
python gltDishMap.py Epr.dat --mask-file mask32.dat --prm-file withphase_aber.prm \
       --output holography_map.pdf --vmin -300 --vmax 300 --cmap seismic
```

**With interpolation options:**
```bash
python gltDishMap.py Epr.dat --holo-smooth bilinear --holo-filterrad 0.8 \
       --output holography_map.pdf --vmin -300 --vmax 300 --cmap seismic
```

### Verification Steps

1. **Check pixel spacing**:
   - Output should show: `Pixel spacing: 124.625 mm`
   - Field of view: `15.952 m`

2. **Check mask**:
   - Output should show: `Using mask from: mask32.dat`
   - Masked pixels: ~3800-4000 of 16384

3. **Visual checks**:
   - Panel boundaries should align with surface features
   - Quadrupod legs should be narrow vertical/horizontal lines
   - Central hole should be small (~750mm diameter)

4. **Compare with FORTRAN**:
   - Panel alignment should match
   - Masked regions should match
   - Color scale values should be comparable

### Files Modified
- `gltDishMap.py` - Updated load_surface_file(), added read_holography_params()

### Files Required
- Input: `Epr.dat` (residual phase from holography)
- Input: `mask32.dat` (mask from holography processing)
- Input: `withphase_aber.prm` (parameters from holography)
- Input: `panelplt.prm` (panel geometry - unchanged)
- Output: PDF map with corrected scale and masking

### Backward Compatibility
- Photogrammetry mode (4-column files) unchanged
- Default parameters provided if .prm file not found
- Falls back to simple circular mask if mask file not available
