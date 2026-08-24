#!/usr/bin/env python3
"""
holoMapGUI.py - GUI interface for GLT Holography Data Reduction Pipeline

A graphical interface that directly runs the holography reduction pipeline,
replacing holoMap.sh with interactive settings for each processing step.

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
from pathlib import Path
from datetime import datetime


class HoloMapGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("GLT Holography Reduction Pipeline")
        self.root.geometry("800x750")
        self.root.minsize(700, 600)
        
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
        
        # === Display Variables ===
        self.no_plot = tk.BooleanVar(value=True)
        self.show_boresight_plot = tk.BooleanVar(value=False)
        
        # === Process tracking ===
        self.process = None
        self.running = False
        self.stop_requested = False
        
        # Build UI
        self.create_widgets()
        
        # Set default data directory
        default_data_dir = os.path.expanduser("~/holo/data")
        if os.path.exists(default_data_dir):
            self.last_dir = default_data_dir
        else:
            self.last_dir = os.getcwd()
    
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
        
        # === Tab 3: Output Log ===
        log_tab = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(log_tab, text="Output Log")
        self.create_log_tab(log_tab)
        
        # === Bottom Button Bar (always visible) ===
        btn_frame = ttk.Frame(self.root, padding="5")
        btn_frame.pack(fill="x", side="bottom")
        
        self.run_btn = ttk.Button(
            btn_frame, text="Run Reduction", command=self.run_pipeline
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
        """Create the main settings tab."""
        parent.columnconfigure(1, weight=1)
        
        row = 0
        
        # === Input File ===
        ttk.Label(parent, text="Input File:", font=('TkDefaultFont', 10, 'bold')).grid(
            row=row, column=0, sticky="w", pady=(0, 5))
        row += 1
        
        file_frame = ttk.Frame(parent)
        file_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 15))
        file_frame.columnconfigure(0, weight=1)
        
        self.file_entry = ttk.Entry(file_frame, textvariable=self.input_file, width=70)
        self.file_entry.grid(row=0, column=0, sticky="ew", padx=(0, 5))
        
        ttk.Button(file_frame, text="Browse...", command=self.browse_file).grid(row=0, column=1)
        
        row += 1
        
        # === Map Size ===
        size_frame = ttk.LabelFrame(parent, text="Map Size", padding="10")
        size_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        sizes = [("32x32", "32"), ("64x64", "64"), ("128x128", "128")]
        for i, (text, value) in enumerate(sizes):
            rb = ttk.Radiobutton(size_frame, text=text, variable=self.map_size, value=value)
            rb.grid(row=0, column=i, padx=15)
        
        row += 1
        
        # === Boresight Calibration ===
        bore_frame = ttk.LabelFrame(parent, text="Boresight Calibration", padding="10")
        bore_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        self.bore_check = ttk.Checkbutton(
            bore_frame, text="Apply boresight phase calibration",
            variable=self.do_boresight_cal
        )
        self.bore_check.grid(row=0, column=0, sticky="w")
        
        self.compare_check = ttk.Checkbutton(
            bore_frame, text="Compare with/without calibration (runs pipeline twice)",
            variable=self.do_compare
        )
        self.compare_check.grid(row=1, column=0, sticky="w", padx=(20, 0), pady=(5, 0))
        
        ttk.Label(bore_frame, text="(Configure detailed options in 'Boresight Cal' tab)",
                  foreground="gray").grid(row=2, column=0, sticky="w", pady=(5, 0))
        
        row += 1
        
        # === Display Options ===
        disp_frame = ttk.LabelFrame(parent, text="Display Options", padding="10")
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
        
        # === Info ===
        info_frame = ttk.LabelFrame(parent, text="Pipeline Steps", padding="10")
        info_frame.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        
        steps_text = """1. detect_raster_start.py  - Detect and trim pre-scan data
2. boresight_cal.py        - Phase drift calibration (if enabled)
3. regrid_holo.py          - Regrid data to regular grid
4. fix_missing_cell.py     - Fix missing grid cells (if boresight cal)
5. preprocess.py           - Geometric corrections
6. holis_aber2.py          - FFT and aberration fitting
7. unwrap2d.py             - 2D phase unwrapping
8. holis_aber2.py --unwrap - Final FFT with unwrapped phase
9. glt_dish_map.py         - Generate surface maps
10. create_summary_pdf.py  - Generate summary PDF"""
        
        ttk.Label(info_frame, text=steps_text, font=('Courier', 9), justify="left").grid(
            row=0, column=0, sticky="w")
    
    def create_boresight_tab(self, parent):
        """Create the boresight calibration settings tab."""
        parent.columnconfigure(1, weight=1)
        
        row = 0
        
        # === Beacon Position ===
        pos_frame = ttk.LabelFrame(parent, text="Beacon Position", padding="10")
        pos_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        ttk.Label(pos_frame, text="Center Azimuth (deg):").grid(row=0, column=0, sticky="w", padx=(0, 10))
        ttk.Entry(pos_frame, textvariable=self.center_az, width=12).grid(row=0, column=1, sticky="w")
        
        ttk.Label(pos_frame, text="Center Elevation (deg):").grid(row=0, column=2, sticky="w", padx=(20, 10))
        ttk.Entry(pos_frame, textvariable=self.center_el, width=12).grid(row=0, column=3, sticky="w")
        
        row += 1
        
        # === Timing Parameters ===
        timing_frame = ttk.LabelFrame(parent, text="Timing Parameters", padding="10")
        timing_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        ttk.Label(timing_frame, text="Slew Margin (s):").grid(row=0, column=0, sticky="w", padx=(0, 10))
        slew_entry = ttk.Entry(timing_frame, textvariable=self.slew_margin, width=8)
        slew_entry.grid(row=0, column=1, sticky="w")
        ttk.Label(timing_frame, text="Time to skip after boresight slew (servo settling)",
                  foreground="gray").grid(row=0, column=2, sticky="w", padx=(10, 0))
        
        row += 1
        
        # === Amplitude Filter ===
        amp_frame = ttk.LabelFrame(parent, text="Amplitude Filter", padding="10")
        amp_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        ttk.Label(amp_frame, text="Minimum Amplitude:").grid(row=0, column=0, sticky="w", padx=(0, 10))
        ttk.Entry(amp_frame, textvariable=self.min_amp, width=8).grid(row=0, column=1, sticky="w")
        ttk.Label(amp_frame, text="Reject boresight visits with amplitude below this",
                  foreground="gray").grid(row=0, column=2, sticky="w", padx=(10, 0))
        
        row += 1
        
        # === Spline Smoothing ===
        smooth_frame = ttk.LabelFrame(parent, text="Spline Smoothing", padding="10")
        smooth_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        ttk.Label(smooth_frame, text="Smoothing Factor:").grid(row=0, column=0, sticky="w", padx=(0, 10))
        ttk.Entry(smooth_frame, textvariable=self.smoothing, width=8).grid(row=0, column=1, sticky="w")
        ttk.Label(smooth_frame, text="Higher = smoother spline (less sensitive to noise)",
                  foreground="gray").grid(row=0, column=2, sticky="w", padx=(10, 0))
        
        row += 1
        
        # === Median Filter ===
        med_frame = ttk.LabelFrame(parent, text="Median Filter (Phase Jump Rejection)", padding="10")
        med_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        ttk.Label(med_frame, text="Window Size:").grid(row=0, column=0, sticky="w", padx=(0, 10))
        ttk.Entry(med_frame, textvariable=self.median_window, width=8).grid(row=0, column=1, sticky="w")
        ttk.Label(med_frame, text="Number of neighboring points for median",
                  foreground="gray").grid(row=0, column=2, sticky="w", padx=(10, 0))
        
        ttk.Label(med_frame, text="Threshold (deg):").grid(row=1, column=0, sticky="w", padx=(0, 10), pady=(5, 0))
        ttk.Entry(med_frame, textvariable=self.median_threshold, width=8).grid(row=1, column=1, sticky="w", pady=(5, 0))
        ttk.Label(med_frame, text="Only filter points deviating more than this from median",
                  foreground="gray").grid(row=1, column=2, sticky="w", padx=(10, 0), pady=(5, 0))
        
        row += 1
        
        # === Recommendations ===
        rec_frame = ttk.LabelFrame(parent, text="Recommendations", padding="10")
        rec_frame.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        rec_text = """* Slew margin 2.0s is recommended (servo needs 1-2s to settle after slew)
* Use min-amp 0.4 if boresight amplitudes are ~0.47-0.48
* Smoothing 128 works well for 128x128 maps
* Median filter helps reject atmospheric phase jumps

Tip: Use "Plot Boresight Phase" button to preview phase drift before running
full reduction. If phase range is <20 deg, calibration is unlikely to help."""
        
        ttk.Label(rec_frame, text=rec_text, justify="left").grid(row=0, column=0, sticky="w")
    
    def create_log_tab(self, parent):
        """Create the output log tab."""
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        
        # Scrolled text widget for output
        self.output_text = scrolledtext.ScrolledText(
            parent, wrap=tk.WORD, font=('Courier', 10)
        )
        self.output_text.grid(row=0, column=0, sticky="nsew")
        
        # Button frame
        btn_frame = ttk.Frame(parent)
        btn_frame.grid(row=1, column=0, sticky="ew", pady=(5, 0))
        
        ttk.Button(btn_frame, text="Clear Log", command=self.clear_log).pack(side="left")
        ttk.Button(btn_frame, text="Save Log...", command=self.save_log).pack(side="left", padx=5)
        
        # Configure text tags
        self.output_text.tag_configure("error", foreground="red")
        self.output_text.tag_configure("success", foreground="green")
        self.output_text.tag_configure("info", foreground="blue")
        self.output_text.tag_configure("step", foreground="purple", font=('Courier', 10, 'bold'))
    
    def browse_file(self):
        """Open file browser to select input file."""
        filename = filedialog.askopenfilename(
            title="Select Holography Data File",
            initialdir=self.last_dir,
            filetypes=[
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        if filename:
            self.input_file.set(filename)
            self.last_dir = os.path.dirname(filename)
    
    def log(self, message, tag=None):
        """Add message to output log."""
        self.output_text.insert(tk.END, message + "\n", tag)
        self.output_text.see(tk.END)
        self.root.update_idletasks()
    
    def clear_log(self):
        """Clear the output log."""
        self.output_text.delete(1.0, tk.END)
    
    def save_log(self):
        """Save log to file."""
        filename = filedialog.asksaveasfilename(
            title="Save Log",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            with open(filename, 'w') as f:
                f.write(self.output_text.get(1.0, tk.END))
            self.log(f"Log saved to {filename}", "info")
    
    def validate_inputs(self):
        """Validate user inputs before running."""
        if not self.input_file.get():
            messagebox.showerror("Error", "Please select an input file")
            return False
        
        if not os.path.exists(self.input_file.get()):
            messagebox.showerror("Error", f"Input file not found:\n{self.input_file.get()}")
            return False
        
        if not os.path.exists(self.work_dir):
            messagebox.showerror("Error", f"Working directory not found:\n{self.work_dir}")
            return False
        
        # Validate numeric inputs
        try:
            float(self.slew_margin.get())
            float(self.smoothing.get())
            int(self.median_window.get())
            float(self.median_threshold.get())
            float(self.min_amp.get())
            float(self.center_az.get())
            float(self.center_el.get())
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid numeric value in settings:\n{e}")
            return False
        
        return True
    
    def run_command(self, cmd, step_name):
        """Run a command and return success status."""
        if self.stop_requested:
            return False
        
        self.root.after(0, lambda: self.log(f"\n{'='*60}", "step"))
        self.root.after(0, lambda: self.log(f"Step: {step_name}", "step"))
        self.root.after(0, lambda: self.log(f"Command: {' '.join(cmd)}"))
        self.root.after(0, lambda: self.log('='*60))
        
        self.root.after(0, lambda: self.progress_var.set(f"Running: {step_name}"))
        
        try:
            # Build shell command with mamba environment activation
            # This ensures all scripts run in the correct environment
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
    
    def build_env_command(self, cmd):
        """
        Build a shell command that activates the mamba environment before running.
        """
        # Convert command list to string
        if isinstance(cmd, list):
            cmd_str = ' '.join(cmd)
        else:
            cmd_str = cmd
        
        # Build the full command with environment activation
        env_cmd = (
            f"source $HOME/.mamba_rc && "
            f"mamba activate nimesh_holo && "
            f"{cmd_str}"
        )
        
        return f"bash -c '{env_cmd}'"
    
    def run_pipeline(self):
        """Run the complete reduction pipeline."""
        if not self.validate_inputs():
            return
        
        if self.running:
            messagebox.showwarning("Warning", "Pipeline is already running")
            return
        
        # Switch to log tab
        self.notebook.select(2)
        
        # Clear log
        self.clear_log()
        
        # Update UI state
        self.running = True
        self.stop_requested = False
        self.run_btn.configure(state="disabled")
        self.stop_btn.configure(state="normal")
        self.status_var.set("Running...")
        
        # Run in background thread
        thread = threading.Thread(target=self.execute_pipeline)
        thread.daemon = True
        thread.start()
    
    def execute_pipeline(self):
        """Execute the pipeline steps (called from background thread)."""
        try:
            input_file = self.input_file.get()
            size = self.map_size.get()
            do_bore = self.do_boresight_cal.get()
            do_compare = self.do_compare.get()
            no_plot = self.no_plot.get()
            
            # Extract prefix from filename
            prefix = Path(input_file).stem
            
            self.root.after(0, lambda: self.log(f"Starting reduction: {prefix}"))
            self.root.after(0, lambda: self.log(f"Map size: {size}x{size}"))
            self.root.after(0, lambda: self.log(f"Boresight cal: {do_bore}"))
            self.root.after(0, lambda: self.log(f"Working directory: {self.work_dir}"))
            
            # Set up parameter file symlinks
            self.setup_symlinks(size)
            
            # === Step 1: Detect raster start ===
            # If interactive mode, we need to handle this specially
            if not no_plot:
                # Run detection in preview mode first
                start_line = self.run_interactive_detection(input_file, do_bore or do_compare)
                if start_line is None:
                    self.root.after(0, lambda: self.log("Detection cancelled by user", "error"))
                    self.root.after(0, self.pipeline_finished)
                    return
                
                # Now run the actual trimming with the confirmed start line
                detect_cmd = ["python", "detect_raster_start.py", input_file, "-o"]
                if do_bore or do_compare:
                    detect_cmd.extend(["trimmed_with_time.txt", "--keep-timestamps"])
                else:
                    detect_cmd.append("trimmed.txt")
                detect_cmd.extend(["--start-line", str(start_line), "--no-plot"])
                
                if not self.run_command(detect_cmd, "Trim data at confirmed start"):
                    self.root.after(0, self.pipeline_finished)
                    return
            else:
                # Non-interactive mode - auto-accept
                detect_cmd = ["python", "detect_raster_start.py", input_file, "-o"]
                if do_bore or do_compare:
                    detect_cmd.extend(["trimmed_with_time.txt", "--keep-timestamps"])
                else:
                    detect_cmd.append("trimmed.txt")
                detect_cmd.append("--no-plot")
                
                if not self.run_command(detect_cmd, "Detect raster start"):
                    self.root.after(0, self.pipeline_finished)
                    return
            
            # === Step 2: Boresight calibration ===
            if do_bore or do_compare:
                bore_cmd = [
                    "python", "boresight_cal.py", "trimmed_with_time.txt",
                    "-o", "calibrated.txt",
                    "--verbose",
                    "--remove-boresight",
                    "--center-az", self.center_az.get(),
                    "--center-el", self.center_el.get(),
                    "--slew-margin", self.slew_margin.get(),
                    "--smoothing", self.smoothing.get(),
                    "--median-window", self.median_window.get(),
                    "--median-threshold", self.median_threshold.get(),
                    "--save-boresight-data",
                    "--save-prefix", prefix
                ]
                
                # Add --min-amp if boresight_cal.py supports it
                try:
                    result = subprocess.run(
                        ["python", "boresight_cal.py", "--help"],
                        capture_output=True, text=True, cwd=self.work_dir
                    )
                    if "--min-amp" in result.stdout:
                        bore_cmd.extend(["--min-amp", self.min_amp.get()])
                except:
                    pass
                
                if self.show_boresight_plot.get() and not no_plot:
                    bore_cmd.append("--plot")
                
                if not do_bore and do_compare:
                    # Only for comparison - don't apply correction
                    bore_cmd.append("--no-correction")
                
                if not self.run_command(bore_cmd, "Boresight calibration"):
                    self.root.after(0, self.pipeline_finished)
                    return
                
                # Extract az/el/amp/phase columns (skip timestamp and comments)
                self.run_command(
                    ["bash", "-c", "awk '!/^#/ {print $2, $3, $4, $5}' calibrated.txt > trimmed.txt"],
                    "Extract calibrated data"
                )
            
            # === Step 3: Regridding ===
            regrid_cmd = ["python", "regrid_holo.py", "trimmed.txt", "regrid.prm"]
            if no_plot:
                regrid_cmd.append("--no-plot")
            
            if not self.run_command(regrid_cmd, "Regrid data"):
                self.root.after(0, self.pipeline_finished)
                return
            
            # === Step 4: Fix missing cells ===
            if do_bore or do_compare:
                if not self.run_command(
                    ["python", "fix_missing_cell.py", size, "--force"],
                    "Fix missing grid cells"
                ):
                    self.root.after(0, self.pipeline_finished)
                    return
            
            # === Step 5: Preprocessing ===
            if not self.run_command(["python", "preprocess.py"], "Preprocessing"):
                self.root.after(0, self.pipeline_finished)
                return
            
            # === Step 6: FFT and aberration fitting ===
            if not self.run_command(["python", "holis_aber2.py"], "FFT and aberration fitting"):
                self.root.after(0, self.pipeline_finished)
                return
            
            # === Step 7: Phase unwrapping ===
            unwrap_cmd = ["python", "unwrap2d.py", "-d", size]
            if not no_plot:
                unwrap_cmd.append("--plot")
            if not self.run_command(unwrap_cmd, "Phase unwrapping"):
                self.root.after(0, self.pipeline_finished)
                return
            
            # === Step 8: Final FFT ===
            if not self.run_command(["python", "holis_aber2.py", "--unwrap"], "Final FFT"):
                self.root.after(0, self.pipeline_finished)
                return
            
            # === Step 9: Copy results ===
            self.root.after(0, lambda: self.log("\nCopying results to results/ directory...", "info"))
            os.makedirs(os.path.join(self.work_dir, "results"), exist_ok=True)
            
            copy_cmds = [
                f"cp ampout.dat results/{prefix}.ampout",
                f"cp phaseout.dat results/{prefix}.phaseout",
                f"cp rgin.dat results/{prefix}.rgrd",
                f"cp Epr.dat results/{prefix}_Epr.dat",
                f"cp holis.log results/{prefix}.log",
                f"cp Ep_um.dat results/{prefix}.Ep_um.dat 2>/dev/null || true",
                f"cp Ea_um.dat results/{prefix}.Ea_um.dat 2>/dev/null || true",
            ]
            
            for cmd in copy_cmds:
                subprocess.run(cmd, shell=True, cwd=self.work_dir)
            
            # === Step 10: Generate maps ===
            # Check which mask file exists
            mask_file = f"mask{size}.dat"
            mask_path = os.path.join(self.work_dir, mask_file)
            if not os.path.exists(mask_path):
                # Try mask32.dat as fallback
                if os.path.exists(os.path.join(self.work_dir, "mask32.dat")):
                    mask_file = "mask32.dat"
                    self.root.after(0, lambda: self.log(f"Note: mask{size}.dat not found, using mask32.dat", "info"))
                else:
                    self.root.after(0, lambda: self.log(f"Warning: No mask file found, skipping mask", "error"))
                    mask_file = None
            
            # Illumination map
            ea_file = os.path.join(self.work_dir, f"results/{prefix}.Ea_um.dat")
            if os.path.exists(ea_file):
                # Show interactively first (if not no_plot)
                if not no_plot:
                    self.run_command(
                        ["python", "glt_dish_map.py", f"results/{prefix}.Ea_um.dat",
                         "--x-shift", "-65", "--y-shift", "65"],
                        "Display illumination map"
                    )
                # Then save to PDF
                self.run_command(
                    ["python", "glt_dish_map.py", f"results/{prefix}.Ea_um.dat",
                     "--x-shift", "-65", "--y-shift", "65",
                     "--output", f"results/{prefix}_illumination.pdf"],
                    "Save illumination map to PDF"
                )
            
            # Surface error map - show interactively first (if not no_plot)
            if not no_plot:
                map_cmd_display = [
                    "python", "glt_dish_map.py", f"results/{prefix}_Epr.dat",
                    "--vmin", "-180", "--vmax", "180",
                    "--x-shift", "-65", "--y-shift", "65",
                    "--prm-file", "withphase_aber.prm"
                ]
                if mask_file:
                    map_cmd_display.extend(["--mask-file", mask_file])
                self.run_command(map_cmd_display, "Display surface error map")
            
            # Then save to PDF
            map_cmd = [
                "python", "glt_dish_map.py", f"results/{prefix}_Epr.dat",
                "--vmin", "-180", "--vmax", "180",
                "--x-shift", "-65", "--y-shift", "65",
                "--prm-file", "withphase_aber.prm",
                "--output", f"results/{prefix}_surface_map.pdf"
            ]
            if mask_file:
                map_cmd.extend(["--mask-file", mask_file])
            
            if not self.run_command(map_cmd, "Save surface map to PDF"):
                self.root.after(0, self.pipeline_finished)
                return
            
            # === Step 11: Summary PDF ===
            # Ask for optional comment (if not in no_plot mode)
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
                comment_event.wait(timeout=300)  # 5 minute timeout
                user_comment = comment_holder['value']
            
            # Build summary command
            summary_cmd = ["python", "create_summary_pdf.py", input_file, f"results/{prefix}"]
            if user_comment:
                summary_cmd.extend(["--comment", user_comment])
            
            if not self.run_command(summary_cmd, "Generate summary PDF"):
                self.root.after(0, lambda: self.log("Warning: Summary PDF generation failed", "error"))
            
            # === Comparison mode ===
            if do_compare:
                self.root.after(0, lambda: self.log("\n" + "="*60, "step"))
                self.root.after(0, lambda: self.log("COMPARISON: Running WITHOUT calibration", "step"))
                self.root.after(0, lambda: self.log("="*60, "step"))
                
                # Save calibrated results
                subprocess.run(f"cp Epr.dat results/{prefix}_cal_Epr.dat", shell=True, cwd=self.work_dir)
                
                # Run without calibration
                bore_cmd_nocal = [
                    "python", "boresight_cal.py", "trimmed_with_time.txt",
                    "-o", "calibrated_nocal.txt",
                    "--verbose", "--no-correction",
                    "--slew-margin", self.slew_margin.get()
                ]
                
                if not self.run_command(bore_cmd_nocal, "Remove boresight (no cal)"):
                    self.root.after(0, self.pipeline_finished)
                    return
                
                self.run_command(
                    ["bash", "-c", "awk '!/^#/ {print $2, $3, $4, $5}' calibrated_nocal.txt > trimmed.txt"],
                    "Extract uncalibrated data"
                )
                
                # Re-run pipeline steps
                for step_name, cmd in [
                    ("Regrid (nocal)", ["python", "regrid_holo.py", "trimmed.txt", "regrid.prm", "--no-plot"]),
                    ("Fix missing cells (nocal)", ["python", "fix_missing_cell.py", size, "--force"]),
                    ("Preprocessing (nocal)", ["python", "preprocess.py"]),
                    ("FFT (nocal)", ["python", "holis_aber2.py"]),
                    ("Phase unwrap (nocal)", ["python", "unwrap2d.py", "-d", size]),
                    ("Final FFT (nocal)", ["python", "holis_aber2.py", "--unwrap"]),
                ]:
                    if not self.run_command(cmd, step_name):
                        self.root.after(0, self.pipeline_finished)
                        return
                
                # Save uncalibrated results
                subprocess.run(f"cp Epr.dat results/{prefix}_nocal_Epr.dat", shell=True, cwd=self.work_dir)
                
                # Generate comparison plot
                self.generate_comparison_plot(prefix)
            
            # === Done ===
            self.root.after(0, lambda: self.log("\n" + "="*60, "success"))
            self.root.after(0, lambda: self.log("REDUCTION COMPLETED SUCCESSFULLY", "success"))
            self.root.after(0, lambda: self.log("="*60, "success"))
            self.root.after(0, lambda: self.log(f"\nResults saved to: {self.work_dir}/results/"))
            
        except Exception as e:
            self.root.after(0, lambda: self.log(f"\nError: {str(e)}", "error"))
            import traceback
            self.root.after(0, lambda: self.log(traceback.format_exc(), "error"))
        
        finally:
            self.root.after(0, self.pipeline_finished)
    
    def setup_symlinks(self, size):
        """Set up parameter file symlinks."""
        symlinks = [
            (f"regrid_{size}x{size}.prm", "regrid.prm"),
            (f"preprocess_{size}.prm", "preprocess.prm"),
            (f"withphase_aber_{size}.prm", "withphase_aber.prm"),
        ]
        
        for src, dst in symlinks:
            dst_path = os.path.join(self.work_dir, dst)
            src_path = os.path.join(self.work_dir, src)
            
            if os.path.exists(src_path):
                if os.path.islink(dst_path):
                    os.unlink(dst_path)
                elif os.path.exists(dst_path):
                    os.remove(dst_path)
                os.symlink(src, dst_path)
    
    def generate_comparison_plot(self, prefix):
        """Generate comparison plot for cal vs nocal."""
        self.root.after(0, lambda: self.log("\nGenerating comparison plot...", "info"))
        
        # Write comparison script to temp file to avoid shell escaping issues
        compare_script = f'''#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

cal_file = "results/{prefix}_cal_Epr.dat"
nocal_file = "results/{prefix}_nocal_Epr.dat"

cal_data = np.loadtxt(cal_file)
nocal_data = np.loadtxt(nocal_file)

grid_size = int(np.sqrt(len(cal_data)))

cal_grid = cal_data.reshape((grid_size, grid_size))
nocal_grid = nocal_data.reshape((grid_size, grid_size))

cal_grid = np.ma.masked_where(cal_grid < -9000, cal_grid)
nocal_grid = np.ma.masked_where(nocal_grid < -9000, nocal_grid)

error_scaling = 299792458.0 / (94.5e9) * 1e6 / (4.0 * np.pi)
cal_grid = cal_grid * error_scaling
nocal_grid = nocal_grid * error_scaling

cal_rms = np.std(cal_grid.compressed())
nocal_rms = np.std(nocal_grid.compressed())

diff_grid = cal_grid - nocal_grid

vmax = max(np.abs(cal_grid).max(), np.abs(nocal_grid).max())

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

im0 = axes[0].imshow(cal_grid, cmap="coolwarm", vmin=-vmax, vmax=vmax, origin="upper")
axes[0].set_title(f"WITH boresight cal\\nRMS: {{cal_rms:.1f}} um")
plt.colorbar(im0, ax=axes[0], label="Surface error (um)")

im1 = axes[1].imshow(nocal_grid, cmap="coolwarm", vmin=-vmax, vmax=vmax, origin="upper")
axes[1].set_title(f"WITHOUT boresight cal\\nRMS: {{nocal_rms:.1f}} um")
plt.colorbar(im1, ax=axes[1], label="Surface error (um)")

diff_max = np.abs(diff_grid).max()
im2 = axes[2].imshow(diff_grid, cmap="coolwarm", vmin=-diff_max, vmax=diff_max, origin="upper")
axes[2].set_title(f"Difference (CAL - NOCAL)\\nMax: {{diff_max:.1f}} um")
plt.colorbar(im2, ax=axes[2], label="Difference (um)")

plt.tight_layout()
plt.savefig(f"results/{prefix}_boresight_comparison.png", dpi=150)
print(f"Saved: results/{prefix}_boresight_comparison.png")
print(f"WITH cal RMS: {{cal_rms:.2f}} um")
print(f"WITHOUT cal RMS: {{nocal_rms:.2f}} um")
plt.show()
'''
        
        try:
            # Write script to temp file
            script_file = os.path.join(self.work_dir, "_temp_compare.py")
            with open(script_file, 'w') as f:
                f.write(compare_script)
            
            # Run with environment
            shell_cmd = self.build_env_command(f"python {script_file}")
            result = subprocess.run(
                shell_cmd,
                capture_output=True, text=True, cwd=self.work_dir, shell=True
            )
            self.root.after(0, lambda: self.log(result.stdout))
            if result.stderr:
                self.root.after(0, lambda: self.log(result.stderr, "error"))
            
            # Clean up temp file
            os.remove(script_file)
        except Exception as e:
            self.root.after(0, lambda: self.log(f"Comparison plot error: {e}", "error"))
    
    def run_interactive_detection(self, input_file, keep_timestamps):
        """
        Run detection in interactive mode with GUI dialogs.
        Returns the confirmed start line, or None if cancelled.
        """
        import tkinter.simpledialog as simpledialog
        
        self.root.after(0, lambda: self.log("\n" + "="*60, "step"))
        self.root.after(0, lambda: self.log("Step: Interactive detection of raster start", "step"))
        self.root.after(0, lambda: self.log("="*60))
        
        # First, run detection to get the auto-detected line
        detect_cmd = f"python detect_raster_start.py {input_file} --detect-only"
        shell_cmd = self.build_env_command(detect_cmd)
        
        try:
            result = subprocess.run(
                shell_cmd,
                capture_output=True, text=True, cwd=self.work_dir, shell=True
            )
            
            self.root.after(0, lambda: self.log(result.stdout))
            
            # Parse detected start line from output
            detected_line = None
            for line in result.stdout.split('\n'):
                if "DETECTED_START:" in line:
                    # Extract the number after the tag
                    import re
                    match = re.search(r'DETECTED_START:\s*(\d+)', line)
                    if match:
                        detected_line = int(match.group(1))
                        break
            
            if detected_line is None:
                self.root.after(0, lambda: self.log("Could not parse detected start line", "error"))
                detected_line = 0
            
            self.root.after(0, lambda: self.log(f"Auto-detected start: line {detected_line}"))
            
            # Show the trajectory plot in separate window FIRST
            plot_cmd = f"python detect_raster_start.py {input_file} --preview-plot {detected_line}"
            plot_shell_cmd = self.build_env_command(plot_cmd)
            
            self.root.after(0, lambda: self.log("Showing trajectory plot..."))
            
            # Run plot in background (non-blocking)
            subprocess.Popen(plot_shell_cmd, cwd=self.work_dir, shell=True)
            
            # Use threading event for synchronization
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
                    # Ask for manual input
                    manual = simpledialog.askinteger(
                        "Enter Start Line",
                        f"Current: {detected_line}\n\nEnter new start line number:",
                        initialvalue=detected_line,
                        minvalue=0
                    )
                    result_holder['value'] = manual  # Will be None if cancelled
                else:
                    result_holder['value'] = None
                
                dialog_event.set()
            
            # Wait 2 seconds for plot to appear, then show dialog
            self.root.after(2000, show_dialog)
            
            # Wait for dialog to complete
            dialog_event.wait(timeout=300)  # 5 minute timeout
            
            if result_holder['value'] == 'WAITING':
                return None  # Timeout
            
            return result_holder['value']
            
        except Exception as e:
            self.root.after(0, lambda: self.log(f"Detection error: {str(e)}", "error"))
            return None
    
    def pipeline_finished(self):
        """Called when pipeline finishes."""
        self.running = False
        self.process = None
        self.run_btn.configure(state="normal")
        self.stop_btn.configure(state="disabled")
        self.status_var.set("Ready")
        self.progress_var.set("")
    
    def stop_pipeline(self):
        """Stop the running pipeline."""
        self.stop_requested = True
        if self.process:
            self.process.terminate()
        self.log("\nPipeline stopped by user", "error")
    
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
        
        self.notebook.select(2)  # Switch to log tab
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
    
    # Try to use a modern theme
    try:
        style = ttk.Style()
        if 'clam' in style.theme_names():
            style.theme_use('clam')
    except:
        pass
    
    app = HoloMapGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
