#!/usr/bin/env python3
"""
holoMapGUI.py - GUI interface for GLT Holography Data Reduction Pipeline

A graphical interface that directly runs the holography reduction pipeline,
replacing holoMap.sh with interactive settings for each processing step.

Features:
- Parameter editing for preprocess.prm and withphase_aber.prm
- Iterative workflow: re-run from any step without starting over
- Individual step buttons for exploratory parameter tuning

Usage:
    python holoMapGUI.py

Requires: tkinter (built into Python)
"""

import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import subprocess
import threading
import os
import sys
import re
from pathlib import Path
from datetime import datetime


class HoloMapGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("GLT Holography Reduction Pipeline")
        self.root.geometry("900x800")
        self.root.minsize(800, 700)
        
        # Working directory
        self.work_dir = os.path.expanduser("~/holo/holoReducePy")
        
        # === Basic Variables ===
        self.input_file = tk.StringVar()
        self.map_size = tk.StringVar(value="128")
        
        # === Boresight Cal Variables ===
        self.do_boresight_cal = tk.BooleanVar(value=False)
        self.do_compare = tk.BooleanVar(value=False)
        self.slew_margin = tk.StringVar(value="2.0")
        self.smoothing = tk.StringVar(value="128")
        self.median_window = tk.StringVar(value="5")
        self.median_threshold = tk.StringVar(value="10")
        self.min_amp = tk.StringVar(value="0.5")
        self.center_az = tk.StringVar(value="230.179")
        self.center_el = tk.StringVar(value="2.8724")
        
        # === Preprocess Parameters ===
        self.ref_plane_dist = tk.StringVar(value="7.48")
        self.do_taper = tk.BooleanVar(value=False)
        self.do_interp = tk.BooleanVar(value=True)
        
        # === Withphase Aber Parameters ===
        self.tx_range = tk.StringVar(value="2820.0")
        self.defocus_corr = tk.StringVar(value="0.0096")
        self.fit_dc = tk.BooleanVar(value=True)
        self.fit_tilt_x = tk.BooleanVar(value=True)
        self.fit_tilt_y = tk.BooleanVar(value=True)
        self.fit_defocus = tk.BooleanVar(value=True)
        self.fit_astig45 = tk.BooleanVar(value=False)
        self.fit_astig = tk.BooleanVar(value=False)
        self.fit_coma_x = tk.BooleanVar(value=True)
        self.fit_coma_y = tk.BooleanVar(value=True)
        
        # === Display Variables ===
        self.no_plot = tk.BooleanVar(value=True)
        self.show_boresight_plot = tk.BooleanVar(value=False)
        
        # === Process tracking ===
        self.process = None
        self.running = False
        self.stop_requested = False
        
        # === Step completion tracking ===
        self.steps_completed = {
            'detect': False,
            'boresight': False,
            'regrid': False,
            'preprocess': False,
            'fft1': False,
            'unwrap': False,
            'fft2': False,
            'visualize': False
        }
        
        # Build UI
        self.create_widgets()
        
        # Set default data directory
        default_data_dir = os.path.expanduser("~/holo/data")
        if os.path.exists(default_data_dir):
            self.last_dir = default_data_dir
        else:
            self.last_dir = os.getcwd()
        
        # Load current parameters from files
        self.load_parameters()
    
    def create_widgets(self):
        """Create all GUI widgets."""
        
        # Create notebook for tabbed interface
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=5, pady=5)
        
        # === Tab 1: Main Settings ===
        main_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(main_tab, text="Main")
        self.create_main_tab(main_tab)
        
        # === Tab 2: Boresight Cal Settings ===
        bore_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(bore_tab, text="Boresight Cal")
        self.create_boresight_tab(bore_tab)
        
        # === Tab 3: Preprocess Parameters ===
        preproc_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(preproc_tab, text="Preprocess")
        self.create_preprocess_tab(preproc_tab)
        
        # === Tab 4: Withphase Aber Parameters ===
        aber_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(aber_tab, text="Aberration Fit")
        self.create_aber_tab(aber_tab)
        
        # === Tab 5: Output Log ===
        log_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(log_tab, text="Output Log")
        self.create_log_tab(log_tab)
        
        # === Bottom Button Bar (always visible) ===
        btn_frame = ttk.Frame(self.root, padding="5")
        btn_frame.pack(fill="x", side="bottom")
        
        self.run_btn = ttk.Button(
            btn_frame, text="Run Full Reduction", command=self.run_pipeline
        )
        self.run_btn.pack(side="left", padx=5)
        
        self.stop_btn = ttk.Button(
            btn_frame, text="Stop", command=self.stop_pipeline, state="disabled"
        )
        self.stop_btn.pack(side="left", padx=5)
        
        ttk.Separator(btn_frame, orient="vertical").pack(side="left", fill="y", padx=10)
        
        self.plot_bore_btn = ttk.Button(
            btn_frame, text="Plot Boresight Phase", command=self.plot_boresight_phase
        )
        self.plot_bore_btn.pack(side="left", padx=5)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        self.progress_var = tk.StringVar(value="")
        
        status_frame = ttk.Frame(self.root)
        status_frame.pack(fill="x", side="bottom")
        
        ttk.Label(status_frame, textvariable=self.status_var, relief="sunken", width=30).pack(side="left", fill="x", expand=True)
        ttk.Label(status_frame, textvariable=self.progress_var, relief="sunken", width=40).pack(side="right")
    
    def create_main_tab(self, parent):
        """Create the main settings tab with step buttons."""
        # Create a canvas with scrollbar for the main tab
        canvas = tk.Canvas(parent)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Bind mouse wheel
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        main_frame = scrollable_frame
        main_frame.columnconfigure(1, weight=1)
        
        row = 0
        
        # === Input File ===
        ttk.Label(main_frame, text="Input File:", font=('TkDefaultFont', 10, 'bold')).grid(
            row=row, column=0, sticky="w", pady=(0, 5))
        row += 1
        
        file_frame = ttk.Frame(main_frame)
        file_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 15))
        file_frame.columnconfigure(0, weight=1)
        
        self.file_entry = ttk.Entry(file_frame, textvariable=self.input_file, width=70)
        self.file_entry.grid(row=0, column=0, sticky="ew", padx=(0, 5))
        
        ttk.Button(file_frame, text="Browse...", command=self.browse_file).grid(row=0, column=1)
        
        row += 1
        
        # === Map Size ===
        size_frame = ttk.LabelFrame(main_frame, text="Map Size", padding="10")
        size_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        sizes = [("32x32", "32"), ("64x64", "64"), ("128x128", "128")]
        for i, (text, value) in enumerate(sizes):
            rb = ttk.Radiobutton(size_frame, text=text, variable=self.map_size, value=value)
            rb.grid(row=0, column=i, padx=15)
        
        row += 1
        
        # === Boresight Calibration Quick Settings ===
        bore_frame = ttk.LabelFrame(main_frame, text="Boresight Calibration", padding="10")
        bore_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        self.bore_check = ttk.Checkbutton(
            bore_frame, text="Apply boresight phase calibration",
            variable=self.do_boresight_cal
        )
        self.bore_check.grid(row=0, column=0, sticky="w")
        
        self.compare_check = ttk.Checkbutton(
            bore_frame, text="Compare with/without calibration",
            variable=self.do_compare
        )
        self.compare_check.grid(row=1, column=0, sticky="w", padx=(20, 0), pady=(5, 0))
        
        ttk.Label(bore_frame, text="(Configure detailed options in 'Boresight Cal' tab)",
                  foreground="gray").grid(row=2, column=0, sticky="w", pady=(5, 0))
        
        row += 1
        
        # === Display Options ===
        disp_frame = ttk.LabelFrame(main_frame, text="Display Options", padding="10")
        disp_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        ttk.Checkbutton(
            disp_frame, text="Skip interactive plots (faster, recommended for remote)",
            variable=self.no_plot
        ).grid(row=0, column=0, sticky="w")
        
        ttk.Checkbutton(
            disp_frame, text="Show boresight calibration diagnostic plot",
            variable=self.show_boresight_plot
        ).grid(row=1, column=0, sticky="w", pady=(5, 0))
        
        row += 1
        
        # === Individual Step Buttons ===
        steps_frame = ttk.LabelFrame(main_frame, text="Pipeline Steps (Run Individually)", padding="10")
        steps_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        # Create step buttons with status indicators
        self.step_buttons = {}
        self.step_status = {}
        
        steps = [
            ("1. Detect & Trim", "detect", self.run_step_detect),
            ("2. Boresight Cal", "boresight", self.run_step_boresight),
            ("3. Regrid", "regrid", self.run_step_regrid),
            ("4. Preprocess", "preprocess", self.run_step_preprocess),
            ("5. FFT + Aberration", "fft1", self.run_step_fft1),
            ("6. Phase Unwrap", "unwrap", self.run_step_unwrap),
            ("7. FFT (unwrapped)", "fft2", self.run_step_fft2),
            ("8. Visualize & Save", "visualize", self.run_step_visualize),
        ]
        
        for i, (label, key, cmd) in enumerate(steps):
            btn_frame_inner = ttk.Frame(steps_frame)
            btn_frame_inner.grid(row=i//2, column=(i%2)*2, columnspan=2, sticky="w", padx=5, pady=3)
            
            # Status indicator (green circle when complete)
            self.step_status[key] = ttk.Label(btn_frame_inner, text="○", width=2)
            self.step_status[key].grid(row=0, column=0)
            
            btn = ttk.Button(btn_frame_inner, text=label, command=cmd, width=18)
            btn.grid(row=0, column=1, padx=(2, 10))
            self.step_buttons[key] = btn
        
        # Add "Run From Here" dropdown
        ttk.Label(steps_frame, text="─" * 60).grid(row=4, column=0, columnspan=4, pady=5)
        
        run_from_frame = ttk.Frame(steps_frame)
        run_from_frame.grid(row=5, column=0, columnspan=4, pady=5)
        
        ttk.Label(run_from_frame, text="Run from step:").grid(row=0, column=0, padx=5)
        
        self.run_from_step = tk.StringVar(value="1. Detect & Trim")
        step_names = [s[0] for s in steps]
        self.run_from_combo = ttk.Combobox(run_from_frame, textvariable=self.run_from_step, 
                                            values=step_names, width=20, state="readonly")
        self.run_from_combo.grid(row=0, column=1, padx=5)
        
        ttk.Button(run_from_frame, text="Run From Selected", 
                   command=self.run_from_selected_step).grid(row=0, column=2, padx=5)
        
        row += 1
        
        # === Quick Info ===
        info_frame = ttk.LabelFrame(main_frame, text="Notes", padding="10")
        info_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        notes_text = """• "Run Full Reduction" runs all steps sequentially (like holoMap.sh)
• Individual step buttons allow re-running from any point
• Green ● indicates step completed in this session
• Modify parameters in Preprocess/Aberration Fit tabs, then re-run that step"""
        
        ttk.Label(info_frame, text=notes_text, justify="left").grid(row=0, column=0, sticky="w")
    
    def create_boresight_tab(self, parent):
        """Create the boresight calibration settings tab."""
        parent.columnconfigure(1, weight=1)
        
        row = 0
        
        # Beacon Position
        pos_frame = ttk.LabelFrame(parent, text="Beacon Position", padding="10")
        pos_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        pos_frame.columnconfigure(1, weight=1)
        pos_frame.columnconfigure(3, weight=1)
        
        ttk.Label(pos_frame, text="Center Az:").grid(row=0, column=0, sticky="e", padx=5)
        ttk.Entry(pos_frame, textvariable=self.center_az, width=12).grid(row=0, column=1, sticky="w")
        ttk.Label(pos_frame, text="deg").grid(row=0, column=2, sticky="w", padx=(2, 20))
        
        ttk.Label(pos_frame, text="Center El:").grid(row=0, column=3, sticky="e", padx=5)
        ttk.Entry(pos_frame, textvariable=self.center_el, width=12).grid(row=0, column=4, sticky="w")
        ttk.Label(pos_frame, text="deg").grid(row=0, column=5, sticky="w")
        
        row += 1
        
        # Calibration Parameters
        cal_frame = ttk.LabelFrame(parent, text="Calibration Parameters", padding="10")
        cal_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        params = [
            ("Slew margin:", self.slew_margin, "seconds", "Time to exclude around boresight"),
            ("Min amplitude:", self.min_amp, "", "Minimum valid amplitude (0-1)"),
            ("Smoothing factor:", self.smoothing, "", "Spline smoothing (higher = smoother)"),
            ("Median window:", self.median_window, "points", "Window for outlier rejection"),
            ("Median threshold:", self.median_threshold, "deg", "Phase jump threshold"),
        ]
        
        for i, (label, var, unit, tooltip) in enumerate(params):
            ttk.Label(cal_frame, text=label).grid(row=i, column=0, sticky="e", padx=5, pady=2)
            e = ttk.Entry(cal_frame, textvariable=var, width=10)
            e.grid(row=i, column=1, sticky="w", pady=2)
            ttk.Label(cal_frame, text=unit).grid(row=i, column=2, sticky="w", padx=5)
            ttk.Label(cal_frame, text=f"({tooltip})", foreground="gray").grid(row=i, column=3, sticky="w")
        
        row += 1
        
        # Recommendations
        rec_frame = ttk.LabelFrame(parent, text="Recommendations", padding="10")
        rec_frame.grid(row=row, column=0, columnspan=2, sticky="nsew", pady=(0, 10))
        parent.rowconfigure(row, weight=1)
        
        rec_text = """When to use boresight calibration:
• Phase range > 50°: Calibration recommended
• Phase range 20-50°: May provide modest improvement
• Phase range < 20°: Calibration unlikely to help

Use 'Plot Boresight Phase' button to analyze phase drift before deciding."""
        
        self.rec_label = ttk.Label(rec_frame, text=rec_text, justify="left")
        self.rec_label.grid(row=0, column=0, sticky="nw")
    
    def create_preprocess_tab(self, parent):
        """Create the preprocess.prm parameters tab."""
        parent.columnconfigure(1, weight=1)
        
        # Header
        ttk.Label(parent, text="Preprocess Parameters (preprocess.prm)", 
                  font=('TkDefaultFont', 11, 'bold')).grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 15))
        
        # Parameters frame
        param_frame = ttk.LabelFrame(parent, text="Geometric Corrections", padding="15")
        param_frame.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        # Reference plane distance
        ttk.Label(param_frame, text="Reference plane distance:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        ttk.Entry(param_frame, textvariable=self.ref_plane_dist, width=12).grid(row=0, column=1, sticky="w", pady=5)
        ttk.Label(param_frame, text="m").grid(row=0, column=2, sticky="w", padx=5)
        ttk.Label(param_frame, text="(Distance from dish vertex to reference horn)", 
                  foreground="gray").grid(row=0, column=3, sticky="w", padx=10)
        
        # Taper option
        ttk.Checkbutton(param_frame, text="Apply far-field taper correction",
                        variable=self.do_taper).grid(row=1, column=0, columnspan=4, sticky="w", pady=5)
        
        # Interpolation option
        ttk.Checkbutton(param_frame, text="Apply interpolation",
                        variable=self.do_interp).grid(row=2, column=0, columnspan=4, sticky="w", pady=5)
        
        # Action buttons
        btn_frame = ttk.Frame(parent)
        btn_frame.grid(row=2, column=0, columnspan=3, pady=15)
        
        ttk.Button(btn_frame, text="Save to preprocess.prm", 
                   command=self.save_preprocess_params).grid(row=0, column=0, padx=5)
        ttk.Button(btn_frame, text="Reload from file", 
                   command=self.load_preprocess_params).grid(row=0, column=1, padx=5)
        ttk.Button(btn_frame, text="Run Preprocess Step", 
                   command=self.run_step_preprocess).grid(row=0, column=2, padx=5)
        
        # Info
        info_frame = ttk.LabelFrame(parent, text="Notes", padding="10")
        info_frame.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(10, 0))
        
        info_text = """After changing parameters:
1. Click "Save to preprocess.prm" to update the parameter file
2. Click "Run Preprocess Step" to re-run from this step
3. The pipeline will continue with FFT and subsequent steps

Input: rgin.dat (from regrid step)
Output: ampout.dat, phaseout.dat"""
        
        ttk.Label(info_frame, text=info_text, justify="left").grid(row=0, column=0, sticky="w")
    
    def create_aber_tab(self, parent):
        """Create the withphase_aber.prm parameters tab."""
        parent.columnconfigure(1, weight=1)
        
        # Header
        ttk.Label(parent, text="Aberration Fitting Parameters (withphase_aber.prm)", 
                  font=('TkDefaultFont', 11, 'bold')).grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 15))
        
        # Physical parameters
        phys_frame = ttk.LabelFrame(parent, text="Physical Parameters", padding="15")
        phys_frame.grid(row=1, column=0, columnspan=4, sticky="ew", pady=(0, 10))
        
        ttk.Label(phys_frame, text="Transmitter range:").grid(row=0, column=0, sticky="e", padx=5, pady=5)
        ttk.Entry(phys_frame, textvariable=self.tx_range, width=12).grid(row=0, column=1, sticky="w", pady=5)
        ttk.Label(phys_frame, text="m").grid(row=0, column=2, sticky="w", padx=5)
        ttk.Label(phys_frame, text="(Distance to beacon)", foreground="gray").grid(row=0, column=3, sticky="w")
        
        ttk.Label(phys_frame, text="Defocus correction:").grid(row=1, column=0, sticky="e", padx=5, pady=5)
        ttk.Entry(phys_frame, textvariable=self.defocus_corr, width=12).grid(row=1, column=1, sticky="w", pady=5)
        ttk.Label(phys_frame, text="mm").grid(row=1, column=2, sticky="w", padx=5)
        ttk.Label(phys_frame, text="(Subreflector defocus offset)", foreground="gray").grid(row=1, column=3, sticky="w")
        
        # Fitting functions
        fit_frame = ttk.LabelFrame(parent, text="Fitting Functions (toggle on/off)", padding="15")
        fit_frame.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(0, 10))
        
        # Create toggle buttons in a grid
        fit_options = [
            ("DC offset", self.fit_dc, "Constant offset"),
            ("Tilt X", self.fit_tilt_x, "Pointing error in X"),
            ("Tilt Y", self.fit_tilt_y, "Pointing error in Y"),
            ("Defocus", self.fit_defocus, "Focus offset"),
            ("Astigmatism-45°", self.fit_astig45, "45° astigmatism"),
            ("Astigmatism", self.fit_astig, "0°/90° astigmatism"),
            ("Coma X", self.fit_coma_x, "Coma in X direction"),
            ("Coma Y", self.fit_coma_y, "Coma in Y direction"),
        ]
        
        for i, (label, var, tooltip) in enumerate(fit_options):
            row_idx = i // 2
            col_idx = (i % 2) * 2
            
            cb = ttk.Checkbutton(fit_frame, text=label, variable=var)
            cb.grid(row=row_idx, column=col_idx, sticky="w", padx=10, pady=3)
            
            ttk.Label(fit_frame, text=f"({tooltip})", foreground="gray").grid(
                row=row_idx, column=col_idx+1, sticky="w", padx=(0, 20))
        
        # Action buttons
        btn_frame = ttk.Frame(parent)
        btn_frame.grid(row=3, column=0, columnspan=4, pady=15)
        
        ttk.Button(btn_frame, text="Save to withphase_aber.prm", 
                   command=self.save_aber_params).grid(row=0, column=0, padx=5)
        ttk.Button(btn_frame, text="Reload from file", 
                   command=self.load_aber_params).grid(row=0, column=1, padx=5)
        ttk.Button(btn_frame, text="Run FFT + Fit Step", 
                   command=self.run_step_fft1).grid(row=0, column=2, padx=5)
        
        # Info
        info_frame = ttk.LabelFrame(parent, text="Notes", padding="10")
        info_frame.grid(row=4, column=0, columnspan=4, sticky="ew", pady=(10, 0))
        
        info_text = """After changing parameters:
1. Click "Save to withphase_aber.prm" to update the parameter file
2. Click "Run FFT + Fit Step" to re-run holis_aber2.py
3. Then run Phase Unwrap and FFT (unwrapped) steps to see final result

Typically enabled: DC, Tilt X/Y, Defocus, Coma X/Y
Typically disabled: Astigmatism terms (unless dish has significant astigmatism)"""
        
        ttk.Label(info_frame, text=info_text, justify="left").grid(row=0, column=0, sticky="w")
    
    def create_log_tab(self, parent):
        """Create the output log tab."""
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        
        # Log text area
        self.log_text = scrolledtext.ScrolledText(
            parent, wrap=tk.WORD, height=25, font=('Courier', 9)
        )
        self.log_text.grid(row=0, column=0, sticky="nsew")
        
        # Configure tags for colored output
        self.log_text.tag_configure("error", foreground="red")
        self.log_text.tag_configure("info", foreground="blue")
        self.log_text.tag_configure("step", foreground="green", font=('Courier', 9, 'bold'))
        
        # Button row
        btn_frame = ttk.Frame(parent)
        btn_frame.grid(row=1, column=0, sticky="ew", pady=(5, 0))
        
        ttk.Button(btn_frame, text="Clear Log", command=self.clear_log).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Save Log...", command=self.save_log).pack(side="left", padx=5)
    
    # ========== Parameter File I/O ==========
    
    def load_parameters(self):
        """Load parameters from both .prm files."""
        self.load_preprocess_params()
        self.load_aber_params()
    
    def load_preprocess_params(self):
        """Load parameters from preprocess.prm."""
        size = self.map_size.get()
        prm_file = os.path.join(self.work_dir, f"preprocess_{size}.prm")
        
        if not os.path.exists(prm_file):
            prm_file = os.path.join(self.work_dir, "preprocess.prm")
        
        if not os.path.exists(prm_file):
            return
        
        try:
            with open(prm_file, 'r') as f:
                content = f.read()
            
            # Parse reference plane distance (dist2)
            match = re.search(r'distance.*reference horn.*?=\s*([\d.]+)', content, re.IGNORECASE)
            if match:
                self.ref_plane_dist.set(match.group(1))
            
            # Parse taper option
            match = re.search(r'do far-field taper.*?=\s*(\d)', content, re.IGNORECASE)
            if match:
                self.do_taper.set(match.group(1) == '1')
            
            # Parse interpolation option
            match = re.search(r'do interpolation.*?=\s*(\d)', content, re.IGNORECASE)
            if match:
                self.do_interp.set(match.group(1) == '1')
                
        except Exception as e:
            self.log(f"Warning: Could not load preprocess params: {e}", "error")
    
    def save_preprocess_params(self):
        """Save parameters to preprocess.prm."""
        size = self.map_size.get()
        prm_file = os.path.join(self.work_dir, f"preprocess_{size}.prm")
        
        if not os.path.exists(prm_file):
            prm_file = os.path.join(self.work_dir, "preprocess.prm")
        
        if not os.path.exists(prm_file):
            messagebox.showerror("Error", f"Parameter file not found: {prm_file}")
            return
        
        try:
            with open(prm_file, 'r') as f:
                lines = f.readlines()
            
            new_lines = []
            for line in lines:
                # Update dist2 (reference horn distance)
                if 'distance' in line.lower() and 'reference horn' in line.lower():
                    # Find the value position (column 50 onwards typically)
                    parts = line.split('=')
                    if len(parts) == 2:
                        prefix = parts[0] + '='
                        line = f"{prefix} {self.ref_plane_dist.get()}\n"
                
                # Update taper option
                elif 'do far-field taper' in line.lower():
                    parts = line.split('=')
                    if len(parts) == 2:
                        prefix = parts[0] + '='
                        line = f"{prefix} {'1' if self.do_taper.get() else '0'}\n"
                
                # Update interpolation option
                elif 'do interpolation' in line.lower():
                    parts = line.split('=')
                    if len(parts) == 2:
                        prefix = parts[0] + '='
                        line = f"{prefix} {'1' if self.do_interp.get() else '0'}\n"
                
                new_lines.append(line)
            
            with open(prm_file, 'w') as f:
                f.writelines(new_lines)
            
            self.log(f"Saved parameters to {prm_file}", "info")
            messagebox.showinfo("Saved", f"Parameters saved to {prm_file}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Could not save parameters: {e}")
    
    def load_aber_params(self):
        """Load parameters from withphase_aber.prm."""
        size = self.map_size.get()
        prm_file = os.path.join(self.work_dir, f"withphase_aber_{size}.prm")
        
        if not os.path.exists(prm_file):
            prm_file = os.path.join(self.work_dir, "withphase_aber.prm")
        
        if not os.path.exists(prm_file):
            return
        
        try:
            with open(prm_file, 'r') as f:
                content = f.read()
            
            # Parse transmitter range
            match = re.search(r'Distance of the transmitter.*?=\s*([\d.]+)', content, re.IGNORECASE)
            if match:
                self.tx_range.set(match.group(1))
            
            # Parse defocus correction
            match = re.search(r'Defocus correction.*?=\s*([\d.]+)', content, re.IGNORECASE)
            if match:
                self.defocus_corr.set(match.group(1))
            
            # Parse fitting functions (1-8)
            fit_vars = [
                self.fit_dc, self.fit_tilt_x, self.fit_tilt_y, self.fit_defocus,
                self.fit_astig45, self.fit_astig, self.fit_coma_x, self.fit_coma_y
            ]
            
            for i, var in enumerate(fit_vars):
                pattern = rf'Fitting function {i+1} ON/OFF.*?=\s*(\d)'
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    var.set(match.group(1) == '1')
                    
        except Exception as e:
            self.log(f"Warning: Could not load aber params: {e}", "error")
    
    def save_aber_params(self):
        """Save parameters to withphase_aber.prm."""
        size = self.map_size.get()
        prm_file = os.path.join(self.work_dir, f"withphase_aber_{size}.prm")
        
        if not os.path.exists(prm_file):
            prm_file = os.path.join(self.work_dir, "withphase_aber.prm")
        
        if not os.path.exists(prm_file):
            messagebox.showerror("Error", f"Parameter file not found: {prm_file}")
            return
        
        try:
            with open(prm_file, 'r') as f:
                lines = f.readlines()
            
            fit_vars = [
                self.fit_dc, self.fit_tilt_x, self.fit_tilt_y, self.fit_defocus,
                self.fit_astig45, self.fit_astig, self.fit_coma_x, self.fit_coma_y
            ]
            
            new_lines = []
            for line in lines:
                # Update transmitter range
                if 'distance of the transmitter' in line.lower():
                    parts = line.split('=')
                    if len(parts) == 2:
                        prefix = parts[0] + '='
                        line = f"{prefix} {self.tx_range.get()}\n"
                
                # Update defocus correction
                elif 'defocus correction' in line.lower():
                    parts = line.split('=')
                    if len(parts) == 2:
                        prefix = parts[0] + '='
                        line = f"{prefix} {self.defocus_corr.get()}\n"
                
                # Update fitting functions
                else:
                    for i, var in enumerate(fit_vars):
                        if f'fitting function {i+1} on/off' in line.lower():
                            parts = line.split('=')
                            if len(parts) == 2:
                                prefix = parts[0] + '='
                                line = f"{prefix} {'1' if var.get() else '0'}\n"
                            break
                
                new_lines.append(line)
            
            with open(prm_file, 'w') as f:
                f.writelines(new_lines)
            
            self.log(f"Saved parameters to {prm_file}", "info")
            messagebox.showinfo("Saved", f"Parameters saved to {prm_file}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Could not save parameters: {e}")
    
    # ========== Logging ==========
    
    def log(self, message, tag=None):
        """Add a message to the log."""
        self.log_text.insert(tk.END, message + "\n", tag)
        self.log_text.see(tk.END)
        self.log_text.update_idletasks()
    
    def clear_log(self):
        """Clear the log text."""
        self.log_text.delete(1.0, tk.END)
    
    def save_log(self):
        """Save log to file."""
        filename = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            initialfile=f"holo_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        )
        if filename:
            with open(filename, 'w') as f:
                f.write(self.log_text.get(1.0, tk.END))
            self.log(f"Log saved to {filename}", "info")
    
    # ========== File Selection ==========
    
    def browse_file(self):
        """Browse for input file."""
        filename = filedialog.askopenfilename(
            initialdir=self.last_dir,
            title="Select holography data file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            self.input_file.set(filename)
            self.last_dir = os.path.dirname(filename)
    
    # ========== Step Status ==========
    
    def update_step_status(self, step, completed):
        """Update the visual status of a step."""
        self.steps_completed[step] = completed
        if step in self.step_status:
            self.step_status[step].config(text="●" if completed else "○",
                                          foreground="green" if completed else "black")
    
    def reset_step_status(self):
        """Reset all step statuses."""
        for key in self.steps_completed:
            self.update_step_status(key, False)
    
    # ========== Command Execution ==========
    
    def build_env_command(self, cmd):
        """Build a shell command that activates the mamba environment."""
        if isinstance(cmd, list):
            cmd_str = ' '.join(cmd)
        else:
            cmd_str = cmd
        
        env_cmd = (
            f"source $HOME/.mamba_rc && "
            f"mamba activate nimesh_holo && "
            f"{cmd_str}"
        )
        
        return f"bash -c '{env_cmd}'"
    
    def run_command(self, cmd, step_name):
        """Run a command and return success status."""
        if self.stop_requested:
            return False
        
        self.root.after(0, lambda: self.log(f"\n{'='*60}", "step"))
        self.root.after(0, lambda: self.log(f"Step: {step_name}", "step"))
        if isinstance(cmd, list):
            self.root.after(0, lambda: self.log(f"Command: {' '.join(cmd)}"))
        else:
            self.root.after(0, lambda: self.log(f"Command: {cmd}"))
        self.root.after(0, lambda: self.log('='*60))
        
        self.root.after(0, lambda: self.progress_var.set(f"Running: {step_name}"))
        
        try:
            shell_cmd = self.build_env_command(cmd)
            
            self.process = subprocess.Popen(
                shell_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=self.work_dir,
                shell=True
            )
            
            for line in self.process.stdout:
                if self.stop_requested:
                    self.process.terminate()
                    return False
                line = line.rstrip()
                self.root.after(0, lambda l=line: self.log(l))
            
            self.process.wait()
            
            if self.process.returncode != 0:
                self.root.after(0, lambda: self.log(f"Error: {step_name} failed with code {self.process.returncode}", "error"))
                return False
            
            return True
            
        except Exception as e:
            self.root.after(0, lambda: self.log(f"Error: {str(e)}", "error"))
            return False
    
    # ========== Individual Step Functions ==========
    
    def run_step_detect(self):
        """Run just the detect/trim step."""
        if not self.input_file.get():
            messagebox.showerror("Error", "Please select an input file")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)  # Switch to log tab
            
            input_file = self.input_file.get()
            no_plot = self.no_plot.get()
            
            # Setup symlinks
            size = self.map_size.get()
            subprocess.run(f"ln -sf regrid_{size}x{size}.prm regrid.prm", shell=True, cwd=self.work_dir)
            subprocess.run(f"ln -sf preprocess_{size}.prm preprocess.prm", shell=True, cwd=self.work_dir)
            subprocess.run(f"ln -sf withphase_aber_{size}.prm withphase_aber.prm", shell=True, cwd=self.work_dir)
            
            if no_plot:
                if self.run_command(
                    f"python detect_raster_start.py {input_file} -o trimmed.txt --no-plot",
                    "Detect and trim raster start"
                ):
                    self.root.after(0, lambda: self.update_step_status('detect', True))
            else:
                # Interactive mode
                start_line = self.run_interactive_detection(input_file, False)
                if start_line is not None:
                    if self.run_command(
                        f"python detect_raster_start.py {input_file} -o trimmed.txt --start-line {start_line} --no-plot",
                        "Trim data at confirmed start"
                    ):
                        self.root.after(0, lambda: self.update_step_status('detect', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_boresight(self):
        """Run just the boresight calibration step."""
        if not os.path.exists(os.path.join(self.work_dir, "trimmed.txt")):
            messagebox.showerror("Error", "trimmed.txt not found. Run Detect step first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            # Need timestamps for boresight cal
            input_file = self.input_file.get()
            if self.do_boresight_cal.get():
                # Re-run detect with timestamps
                self.run_command(
                    f"python detect_raster_start.py {input_file} -o trimmed_with_time.txt --keep-timestamps --no-plot",
                    "Prepare timestamped data"
                )
                
                bore_cmd = [
                    "python", "boresight_cal.py", "trimmed_with_time.txt",
                    "-o", "calibrated.txt", "--verbose",
                    "--center-az", self.center_az.get(),
                    "--center-el", self.center_el.get(),
                    "--slew-margin", self.slew_margin.get(),
                    "--smoothing", self.smoothing.get(),
                    "--median-window", self.median_window.get(),
                    "--median-threshold", self.median_threshold.get()
                ]
                
                if self.show_boresight_plot.get() and not self.no_plot.get():
                    bore_cmd.append("--plot")
                
                if self.run_command(bore_cmd, "Boresight calibration"):
                    # Extract to trimmed.txt
                    self.run_command(
                        "bash -c \"awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt\"",
                        "Extract calibrated data"
                    )
                    self.root.after(0, lambda: self.update_step_status('boresight', True))
            else:
                self.root.after(0, lambda: self.log("Boresight calibration disabled, skipping", "info"))
                self.root.after(0, lambda: self.update_step_status('boresight', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_regrid(self):
        """Run just the regrid step."""
        if not os.path.exists(os.path.join(self.work_dir, "trimmed.txt")):
            messagebox.showerror("Error", "trimmed.txt not found. Run Detect step first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            no_plot = self.no_plot.get()
            
            regrid_cmd = ["python", "regrid_holo.py", "trimmed.txt", "regrid.prm"]
            if no_plot:
                regrid_cmd.append("--no-plot")
            
            if self.run_command(regrid_cmd, "Regrid data"):
                # Fix missing cells if boresight cal was used
                if self.do_boresight_cal.get():
                    size = self.map_size.get()
                    self.run_command(f"python fix_missing_cell.py {size} --force", "Fix missing cells")
                self.root.after(0, lambda: self.update_step_status('regrid', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_preprocess(self):
        """Run just the preprocess step."""
        if not os.path.exists(os.path.join(self.work_dir, "rgin.dat")):
            messagebox.showerror("Error", "rgin.dat not found. Run Regrid step first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            if self.run_command(["python", "preprocess.py"], "Preprocessing"):
                self.root.after(0, lambda: self.update_step_status('preprocess', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_fft1(self):
        """Run just the FFT and aberration fitting step."""
        if not os.path.exists(os.path.join(self.work_dir, "ampout.dat")):
            messagebox.showerror("Error", "ampout.dat not found. Run Preprocess step first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            if self.run_command(["python", "holis_aber2.py"], "FFT and aberration fitting"):
                self.root.after(0, lambda: self.update_step_status('fft1', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_unwrap(self):
        """Run just the phase unwrapping step."""
        if not os.path.exists(os.path.join(self.work_dir, "Ep.dat")):
            messagebox.showerror("Error", "Ep.dat not found. Run FFT step first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            size = self.map_size.get()
            no_plot = self.no_plot.get()
            
            unwrap_cmd = ["python", "unwrap2d.py", "-d", size]
            if not no_plot:
                unwrap_cmd.append("--plot")
            
            if self.run_command(unwrap_cmd, "Phase unwrapping"):
                self.root.after(0, lambda: self.update_step_status('unwrap', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_fft2(self):
        """Run the final FFT with unwrapped phase."""
        if not os.path.exists(os.path.join(self.work_dir, "tk.dat")):
            messagebox.showerror("Error", "tk.dat not found. Run Unwrap step first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            if self.run_command(["python", "holis_aber2.py", "--unwrap"], "Final FFT with unwrapped phase"):
                self.root.after(0, lambda: self.update_step_status('fft2', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_step_visualize(self):
        """Run visualization and save results."""
        if not os.path.exists(os.path.join(self.work_dir, "Epr.dat")):
            messagebox.showerror("Error", "Epr.dat not found. Run FFT steps first.")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            input_file = self.input_file.get()
            prefix = os.path.splitext(os.path.basename(input_file))[0]
            size = self.map_size.get()
            no_plot = self.no_plot.get()
            
            # Create results directory
            os.makedirs(os.path.join(self.work_dir, "results"), exist_ok=True)
            
            # Copy result files
            copy_cmds = [
                f"cp ampout.dat results/{prefix}.ampout",
                f"cp phaseout.dat results/{prefix}.phaseout",
                f"cp rgin.dat results/{prefix}.rgrd",
                f"cp Epr.dat results/{prefix}_Epr.dat",
                f"cp holis.log results/{prefix}.log 2>/dev/null || true",
                f"cp Ep_um.dat results/{prefix}.Ep_um.dat 2>/dev/null || true",
                f"cp Ea_um.dat results/{prefix}.Ea_um.dat 2>/dev/null || true",
            ]
            
            for cmd in copy_cmds:
                subprocess.run(cmd, shell=True, cwd=self.work_dir)
            
            # Check mask file
            mask_file = f"mask{size}.dat"
            if not os.path.exists(os.path.join(self.work_dir, mask_file)):
                if os.path.exists(os.path.join(self.work_dir, "mask32.dat")):
                    mask_file = "mask32.dat"
                else:
                    mask_file = None
            
            # Illumination map
            ea_file = f"results/{prefix}.Ea_um.dat"
            if os.path.exists(os.path.join(self.work_dir, ea_file)):
                if not no_plot:
                    self.run_command(
                        ["python", "glt_dish_map.py", ea_file, "--x-shift", "-65", "--y-shift", "65"],
                        "Display illumination map"
                    )
                self.run_command(
                    ["python", "glt_dish_map.py", ea_file, "--x-shift", "-65", "--y-shift", "65",
                     "--output", f"results/{prefix}_illumination.pdf"],
                    "Save illumination map"
                )
            
            # Surface error map
            map_cmd = ["python", "glt_dish_map.py", f"results/{prefix}_Epr.dat",
                       "--vmin", "-180", "--vmax", "180",
                       "--x-shift", "-65", "--y-shift", "65",
                       "--prm-file", "withphase_aber.prm"]
            if mask_file:
                map_cmd.extend(["--mask-file", mask_file])
            
            if not no_plot:
                self.run_command(map_cmd, "Display surface error map")
            
            map_cmd.extend(["--output", f"results/{prefix}_surface_map.pdf"])
            self.run_command(map_cmd, "Save surface map")
            
            # Ask for comment
            user_comment = None
            if not no_plot:
                comment_event = threading.Event()
                comment_holder = {'value': None}
                
                def ask_comment():
                    import tkinter.simpledialog as simpledialog
                    comment = simpledialog.askstring(
                        "Summary Comments",
                        "Enter optional comments for summary page:\n(Press Cancel or leave empty to skip)",
                        parent=self.root
                    )
                    comment_holder['value'] = comment
                    comment_event.set()
                
                self.root.after(0, ask_comment)
                comment_event.wait(timeout=300)
                user_comment = comment_holder['value']
            
            # Generate summary
            summary_cmd = ["python", "create_summary_pdf.py", input_file, f"results/{prefix}"]
            if user_comment:
                summary_cmd.extend(["--comment", user_comment])
            
            self.run_command(summary_cmd, "Generate summary PDF")
            self.root.after(0, lambda: self.update_step_status('visualize', True))
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
            self.root.after(0, lambda: self.log("\n✓ Visualization complete!", "info"))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    def run_from_selected_step(self):
        """Run pipeline from selected step onwards."""
        step_map = {
            "1. Detect & Trim": 0,
            "2. Boresight Cal": 1,
            "3. Regrid": 2,
            "4. Preprocess": 3,
            "5. FFT + Aberration": 4,
            "6. Phase Unwrap": 5,
            "7. FFT (unwrapped)": 6,
            "8. Visualize & Save": 7,
        }
        
        selected = self.run_from_step.get()
        start_idx = step_map.get(selected, 0)
        
        if not self.input_file.get() and start_idx == 0:
            messagebox.showerror("Error", "Please select an input file")
            return
        
        def run():
            self.running = True
            self.stop_requested = False
            self.root.after(0, lambda: self.run_btn.config(state="disabled"))
            self.root.after(0, lambda: self.stop_btn.config(state="normal"))
            self.root.after(0, lambda: self.status_var.set("Running..."))
            self.notebook.select(4)
            
            step_funcs = [
                self._exec_detect,
                self._exec_boresight,
                self._exec_regrid,
                self._exec_preprocess,
                self._exec_fft1,
                self._exec_unwrap,
                self._exec_fft2,
                self._exec_visualize,
            ]
            
            for i in range(start_idx, len(step_funcs)):
                if self.stop_requested:
                    break
                if not step_funcs[i]():
                    break
            
            self.running = False
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            self.root.after(0, lambda: self.stop_btn.config(state="disabled"))
            self.root.after(0, lambda: self.status_var.set("Ready"))
            self.root.after(0, lambda: self.progress_var.set(""))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()
    
    # ========== Internal step execution (for run_from_selected) ==========
    
    def _exec_detect(self):
        """Execute detect step inline."""
        input_file = self.input_file.get()
        no_plot = self.no_plot.get()
        size = self.map_size.get()
        
        subprocess.run(f"ln -sf regrid_{size}x{size}.prm regrid.prm", shell=True, cwd=self.work_dir)
        subprocess.run(f"ln -sf preprocess_{size}.prm preprocess.prm", shell=True, cwd=self.work_dir)
        subprocess.run(f"ln -sf withphase_aber_{size}.prm withphase_aber.prm", shell=True, cwd=self.work_dir)
        
        if no_plot:
            result = self.run_command(
                f"python detect_raster_start.py {input_file} -o trimmed.txt --no-plot",
                "Detect and trim raster start"
            )
        else:
            start_line = self.run_interactive_detection(input_file, False)
            if start_line is None:
                return False
            result = self.run_command(
                f"python detect_raster_start.py {input_file} -o trimmed.txt --start-line {start_line} --no-plot",
                "Trim data at confirmed start"
            )
        
        if result:
            self.root.after(0, lambda: self.update_step_status('detect', True))
        return result
    
    def _exec_boresight(self):
        """Execute boresight step inline."""
        input_file = self.input_file.get()
        
        if self.do_boresight_cal.get():
            self.run_command(
                f"python detect_raster_start.py {input_file} -o trimmed_with_time.txt --keep-timestamps --no-plot",
                "Prepare timestamped data"
            )
            
            bore_cmd = [
                "python", "boresight_cal.py", "trimmed_with_time.txt",
                "-o", "calibrated.txt", "--verbose",
                "--center-az", self.center_az.get(),
                "--center-el", self.center_el.get(),
                "--slew-margin", self.slew_margin.get(),
                "--smoothing", self.smoothing.get(),
                "--median-window", self.median_window.get(),
                "--median-threshold", self.median_threshold.get()
            ]
            
            if self.show_boresight_plot.get() and not self.no_plot.get():
                bore_cmd.append("--plot")
            
            if not self.run_command(bore_cmd, "Boresight calibration"):
                return False
            
            self.run_command(
                "bash -c \"awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt\"",
                "Extract calibrated data"
            )
        
        self.root.after(0, lambda: self.update_step_status('boresight', True))
        return True
    
    def _exec_regrid(self):
        """Execute regrid step inline."""
        no_plot = self.no_plot.get()
        
        regrid_cmd = ["python", "regrid_holo.py", "trimmed.txt", "regrid.prm"]
        if no_plot:
            regrid_cmd.append("--no-plot")
        
        if not self.run_command(regrid_cmd, "Regrid data"):
            return False
        
        if self.do_boresight_cal.get():
            size = self.map_size.get()
            self.run_command(f"python fix_missing_cell.py {size} --force", "Fix missing cells")
        
        self.root.after(0, lambda: self.update_step_status('regrid', True))
        return True
    
    def _exec_preprocess(self):
        """Execute preprocess step inline."""
        result = self.run_command(["python", "preprocess.py"], "Preprocessing")
        if result:
            self.root.after(0, lambda: self.update_step_status('preprocess', True))
        return result
    
    def _exec_fft1(self):
        """Execute first FFT step inline."""
        result = self.run_command(["python", "holis_aber2.py"], "FFT and aberration fitting")
        if result:
            self.root.after(0, lambda: self.update_step_status('fft1', True))
        return result
    
    def _exec_unwrap(self):
        """Execute unwrap step inline."""
        size = self.map_size.get()
        no_plot = self.no_plot.get()
        
        unwrap_cmd = ["python", "unwrap2d.py", "-d", size]
        if not no_plot:
            unwrap_cmd.append("--plot")
        
        result = self.run_command(unwrap_cmd, "Phase unwrapping")
        if result:
            self.root.after(0, lambda: self.update_step_status('unwrap', True))
        return result
    
    def _exec_fft2(self):
        """Execute second FFT step inline."""
        result = self.run_command(["python", "holis_aber2.py", "--unwrap"], "Final FFT with unwrapped phase")
        if result:
            self.root.after(0, lambda: self.update_step_status('fft2', True))
        return result
    
    def _exec_visualize(self):
        """Execute visualization step inline."""
        input_file = self.input_file.get()
        prefix = os.path.splitext(os.path.basename(input_file))[0]
        size = self.map_size.get()
        no_plot = self.no_plot.get()
        
        os.makedirs(os.path.join(self.work_dir, "results"), exist_ok=True)
        
        # Copy files
        copy_cmds = [
            f"cp ampout.dat results/{prefix}.ampout",
            f"cp phaseout.dat results/{prefix}.phaseout",
            f"cp rgin.dat results/{prefix}.rgrd",
            f"cp Epr.dat results/{prefix}_Epr.dat",
            f"cp holis.log results/{prefix}.log 2>/dev/null || true",
            f"cp Ep_um.dat results/{prefix}.Ep_um.dat 2>/dev/null || true",
            f"cp Ea_um.dat results/{prefix}.Ea_um.dat 2>/dev/null || true",
        ]
        
        for cmd in copy_cmds:
            subprocess.run(cmd, shell=True, cwd=self.work_dir)
        
        mask_file = f"mask{size}.dat"
        if not os.path.exists(os.path.join(self.work_dir, mask_file)):
            if os.path.exists(os.path.join(self.work_dir, "mask32.dat")):
                mask_file = "mask32.dat"
            else:
                mask_file = None
        
        # Illumination
        ea_file = f"results/{prefix}.Ea_um.dat"
        if os.path.exists(os.path.join(self.work_dir, ea_file)):
            if not no_plot:
                self.run_command(
                    ["python", "glt_dish_map.py", ea_file, "--x-shift", "-65", "--y-shift", "65"],
                    "Display illumination map"
                )
            self.run_command(
                ["python", "glt_dish_map.py", ea_file, "--x-shift", "-65", "--y-shift", "65",
                 "--output", f"results/{prefix}_illumination.pdf"],
                "Save illumination map"
            )
        
        # Surface map
        map_cmd = ["python", "glt_dish_map.py", f"results/{prefix}_Epr.dat",
                   "--vmin", "-180", "--vmax", "180",
                   "--x-shift", "-65", "--y-shift", "65",
                   "--prm-file", "withphase_aber.prm"]
        if mask_file:
            map_cmd.extend(["--mask-file", mask_file])
        
        if not no_plot:
            self.run_command(map_cmd, "Display surface error map")
        
        map_cmd_save = map_cmd + ["--output", f"results/{prefix}_surface_map.pdf"]
        self.run_command(map_cmd_save, "Save surface map")
        
        # Comment dialog
        user_comment = None
        if not no_plot:
            comment_event = threading.Event()
            comment_holder = {'value': None}
            
            def ask_comment():
                import tkinter.simpledialog as simpledialog
                comment = simpledialog.askstring(
                    "Summary Comments",
                    "Enter optional comments for summary page:\n(Press Cancel or leave empty to skip)",
                    parent=self.root
                )
                comment_holder['value'] = comment
                comment_event.set()
            
            self.root.after(0, ask_comment)
            comment_event.wait(timeout=300)
            user_comment = comment_holder['value']
        
        summary_cmd = ["python", "create_summary_pdf.py", input_file, f"results/{prefix}"]
        if user_comment:
            summary_cmd.extend(["--comment", user_comment])
        
        self.run_command(summary_cmd, "Generate summary PDF")
        self.root.after(0, lambda: self.update_step_status('visualize', True))
        self.root.after(0, lambda: self.log("\n✓ Reduction complete!", "info"))
        return True
    
    # ========== Interactive Detection ==========
    
    def run_interactive_detection(self, input_file, keep_timestamps):
        """Run detection in interactive mode with GUI dialogs."""
        import tkinter.simpledialog as simpledialog
        
        self.root.after(0, lambda: self.log("\n" + "="*60, "step"))
        self.root.after(0, lambda: self.log("Step: Interactive detection of raster start", "step"))
        self.root.after(0, lambda: self.log("="*60))
        
        detect_cmd = f"python detect_raster_start.py {input_file} --detect-only"
        shell_cmd = self.build_env_command(detect_cmd)
        
        try:
            result = subprocess.run(
                shell_cmd,
                capture_output=True, text=True, cwd=self.work_dir, shell=True
            )
            
            self.root.after(0, lambda: self.log(result.stdout))
            
            detected_line = None
            for line in result.stdout.split('\n'):
                if "DETECTED_START:" in line:
                    match = re.search(r'DETECTED_START:\s*(\d+)', line)
                    if match:
                        detected_line = int(match.group(1))
                        break
            
            if detected_line is None:
                self.root.after(0, lambda: self.log("Could not parse detected start line", "error"))
                detected_line = 0
            
            self.root.after(0, lambda: self.log(f"Auto-detected start: line {detected_line}"))
            self.root.after(0, lambda: self.log("Showing trajectory plot..."))
            
            plot_cmd = f"python detect_raster_start.py {input_file} --preview-plot {detected_line}"
            plot_shell_cmd = self.build_env_command(plot_cmd)
            
            subprocess.Popen(plot_shell_cmd, cwd=self.work_dir, shell=True)
            
            dialog_event = threading.Event()
            result_holder = {'value': 'WAITING'}
            
            def show_dialog():
                response = messagebox.askyesnocancel(
                    "Confirm Start Line",
                    f"Auto-detected start line: {detected_line}\n\n"
                    f"Check the trajectory plot window.\n"
                    f"The red vertical line shows the detected start.\n\n"
                    f"Accept this start line?\n\n"
                    f"• Yes = Accept line {detected_line}\n"
                    f"• No = Enter a different line number\n"
                    f"• Cancel = Abort reduction"
                )
                
                if response is True:
                    result_holder['value'] = detected_line
                elif response is False:
                    manual = simpledialog.askinteger(
                        "Enter Start Line",
                        f"Current: {detected_line}\n\nEnter new start line number:",
                        initialvalue=detected_line,
                        minvalue=0
                    )
                    result_holder['value'] = manual
                else:
                    result_holder['value'] = None
                
                dialog_event.set()
            
            self.root.after(2000, show_dialog)
            dialog_event.wait(timeout=300)
            
            if result_holder['value'] == 'WAITING':
                return None
            
            return result_holder['value']
            
        except Exception as e:
            self.root.after(0, lambda: self.log(f"Detection error: {str(e)}", "error"))
            return None
    
    # ========== Full Pipeline ==========
    
    def run_pipeline(self):
        """Run the full reduction pipeline."""
        if not self.input_file.get():
            messagebox.showerror("Error", "Please select an input file")
            return
        
        self.reset_step_status()
        self.run_from_step.set("1. Detect & Trim")
        self.run_from_selected_step()
    
    def stop_pipeline(self):
        """Stop the running pipeline."""
        self.stop_requested = True
        if self.process:
            self.process.terminate()
        self.status_var.set("Stopped")
    
    def pipeline_finished(self):
        """Called when pipeline finishes."""
        self.running = False
        self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.status_var.set("Ready")
        self.progress_var.set("")
    
    def plot_boresight_phase(self):
        """Quick action to plot boresight phase."""
        if not self.input_file.get():
            messagebox.showerror("Error", "Please select an input file first")
            return
        
        if not os.path.exists(self.input_file.get()):
            messagebox.showerror("Error", f"Input file not found:\n{self.input_file.get()}")
            return
        
        cmd_str = (
            f"python plot_boresight_phase.py "
            f"{self.input_file.get()} "
            f"--center-az {self.center_az.get()} "
            f"--center-el {self.center_el.get()} "
            f"-v"
        )
        shell_cmd = self.build_env_command(cmd_str)
        
        self.notebook.select(4)
        self.log("=" * 60)
        self.log("Plotting boresight phase...", "info")
        self.log(f"Command: {cmd_str}")
        self.log("=" * 60)
        
        def run():
            try:
                result = subprocess.run(shell_cmd, capture_output=True, text=True, 
                                        cwd=self.work_dir, shell=True)
                self.root.after(0, lambda: self.log(result.stdout))
                if result.stderr:
                    self.root.after(0, lambda: self.log(result.stderr, "error"))
            except Exception as e:
                self.root.after(0, lambda: self.log(f"Error: {str(e)}", "error"))
        
        thread = threading.Thread(target=run)
        thread.daemon = True
        thread.start()


def main():
    root = tk.Tk()
    app = HoloMapGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
