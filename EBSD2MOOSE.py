# Nicolò Grilli
# University of Bristol
# 1 Giugno 2023

# convert EBSD into Euler angles file
# compatible with MOOSE framework or UMAT

import sys
import argparse
from ebsd import EBSD
from abaqus_input_file import AbaqusInputFile
from multiphase import Multiphase
from stereographic_triangle import StereographicTriangle
from odf import ODF
from neper import Neper

def run_conversion(args):
    """Run the EBSD to MOOSE/UMAT conversion with the given arguments."""
    if args.neper:
        data = Neper(args.filename)  # Neper needs mesh file instead of EBSD file
        data.parse_mesh_file()
        data.generate_interface_file()
        data.generate_euler_angles_file()
        return

    data = EBSD(args.filename)

    data.parse_ebsd_file()
    data.generate_2D_Euler_angles_map()

    data.plot_EBSD_map(args.frequency, args.nx_min, args.nx_max, args.ny_min, args.ny_max)

    if args.UMAT:
        data.generate_UMAT_Euler_angles_file(args.filename.rsplit('.', maxsplit=1)[0] + '.txt',args.frequency, args.thickness, args.nx_min, args.nx_max, args.ny_min, args.ny_max)
        abaqus_file = AbaqusInputFile(args.filename.rsplit('.', maxsplit=1)[0] + '.inp', data, args.thickness)
        abaqus_file.write_input_file(args.frequency, args.nx_min, args.nx_max, args.ny_min, args.ny_max)
        if args.aster:
            abaqus_file.inp2med()
    else:
        data.generate_MOOSE_Euler_angles_file(args.filename.rsplit('.', maxsplit=1)[0] + '.txt',args.frequency, args.thickness, args.nx_min, args.nx_max, args.ny_min, args.ny_max)
        if args.multiphase:
            phase = Multiphase(args.filename, data)
            phase.generate_MOOSE_phase_file(args.filename.rsplit('.', maxsplit=1)[0] + '_phase.txt',args.frequency, args.thickness, args.nx_min, args.nx_max, args.ny_min, args.ny_max)
        if args.stereographic_projection:
            stereo = StereographicTriangle(data)
            stereo.generate_stereographic_projection()
            stereo.plot_stereographic_projection()
        if args.odf:
            odf = ODF(data, args.N_grain_sample)
            odf.generate_sample()

def main_gui():
    """Launch the PyQt graphical interface."""
    from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                                 QLabel, QLineEdit, QPushButton, QCheckBox, QSpinBox,
                                 QFileDialog, QMessageBox, QGroupBox, QFormLayout,
                                 QTextEdit, QScrollArea)
    from PyQt6.QtCore import Qt, QThread, pyqtSignal
    from PyQt6.QtGui import QFont

    class ConversionWorker(QThread):
        """Worker thread for running the conversion."""
        progress = pyqtSignal(str)
        finished = pyqtSignal(bool, str)

        def __init__(self, args):
            super().__init__()
            self.args = args

        def run(self):
             self.progress.emit("Starting conversion...")
             run_conversion(self.args)
             self.progress.emit("Conversion completed successfully!")
             self.finished.emit(True, "Conversion completed successfully!")

    class EBSD2MOOSEGUI(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("EBSD2MOOSE - EBSD to MOOSE/UMAT Converter")
            self.setMinimumSize(700, 800)
            self.worker = None
            self.init_ui()

        def init_ui(self):
            central_widget = QWidget()
            self.setCentralWidget(central_widget)
            main_layout = QVBoxLayout(central_widget)

            # Title
            title = QLabel("EBSD2MOOSE Converter")
            title.setFont(QFont("Arial", 16, QFont.Weight.Bold))
            title.setAlignment(Qt.AlignmentFlag.AlignCenter)
            main_layout.addWidget(title)

            # Scroll area for form
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll_content = QWidget()
            scroll_layout = QVBoxLayout(scroll_content)

            # Input File Section
            file_group = QGroupBox("Input File")
            file_layout = QHBoxLayout()
            self.file_edit = QLineEdit()
            self.file_edit.setPlaceholderText("Select EBSD file (ctf or ang)...")
            file_btn = QPushButton("Browse...")
            file_btn.clicked.connect(self.browse_file)
            file_layout.addWidget(self.file_edit)
            file_layout.addWidget(file_btn)
            file_group.setLayout(file_layout)
            scroll_layout.addWidget(file_group)

            # Parameters Section
            params_group = QGroupBox("Parameters")
            params_layout = QFormLayout()

            self.freq_spin = QSpinBox()
            self.freq_spin.setRange(1, 1000)
            self.freq_spin.setValue(1)
            params_layout.addRow("Frequency:", self.freq_spin)

            self.thick_spin = QSpinBox()
            self.thick_spin.setRange(1, 1000)
            self.thick_spin.setValue(1)
            params_layout.addRow("Thickness:", self.thick_spin)

            self.nx_min_spin = QSpinBox()
            self.nx_min_spin.setRange(-1, 100000)
            self.nx_min_spin.setValue(-1)
            self.nx_min_spin.setSpecialValueText("Auto (full range)")
            params_layout.addRow("Min nx:", self.nx_min_spin)

            self.nx_max_spin = QSpinBox()
            self.nx_max_spin.setRange(-1, 100000)
            self.nx_max_spin.setValue(-1)
            self.nx_max_spin.setSpecialValueText("Auto (full range)")
            params_layout.addRow("Max nx:", self.nx_max_spin)

            self.ny_min_spin = QSpinBox()
            self.ny_min_spin.setRange(-1, 100000)
            self.ny_min_spin.setValue(-1)
            self.ny_min_spin.setSpecialValueText("Auto (full range)")
            params_layout.addRow("Min ny:", self.ny_min_spin)

            self.ny_max_spin = QSpinBox()
            self.ny_max_spin.setRange(-1, 100000)
            self.ny_max_spin.setValue(-1)
            self.ny_max_spin.setSpecialValueText("Auto (full range)")
            params_layout.addRow("Max ny:", self.ny_max_spin)

            self.n_grain_spin = QSpinBox()
            self.n_grain_spin.setRange(1, 1000000)
            self.n_grain_spin.setValue(1)
            params_layout.addRow("N Grain Sample (ODF):", self.n_grain_spin)

            params_group.setLayout(params_layout)
            scroll_layout.addWidget(params_group)

            # Options Section
            options_group = QGroupBox("Options")
            options_layout = QVBoxLayout()

            self.umat_check = QCheckBox("Generate UMAT format output")
            self.aster_check = QCheckBox("Generate Aster format (requires UMAT)")
            self.multiphase_check = QCheckBox("Generate multiphase file")
            self.stereo_check = QCheckBox("Generate stereographic projection")
            self.odf_check = QCheckBox("Generate ODF (Orientation Distribution Function)")
            self.neper_check = QCheckBox("Use Neper mesh file (instead of EBSD)")

            options_layout.addWidget(self.umat_check)
            options_layout.addWidget(self.aster_check)
            options_layout.addWidget(self.multiphase_check)
            options_layout.addWidget(self.stereo_check)
            options_layout.addWidget(self.odf_check)
            options_layout.addWidget(self.neper_check)

            options_group.setLayout(options_layout)
            scroll_layout.addWidget(options_group)

            scroll.setWidget(scroll_content)
            main_layout.addWidget(scroll)

            # Log output
            log_group = QGroupBox("Log")
            log_layout = QVBoxLayout()
            self.log_text = QTextEdit()
            self.log_text.setReadOnly(True)
            self.log_text.setMaximumHeight(250)
            log_layout.addWidget(self.log_text)
            log_group.setLayout(log_layout)
            main_layout.addWidget(log_group)

            # Buttons
            btn_layout = QHBoxLayout()
            self.run_btn = QPushButton("Run Conversion")
            self.run_btn.setFont(QFont("Arial", 12, QFont.Weight.Bold))
            self.run_btn.setMinimumHeight(40)
            self.run_btn.clicked.connect(self.run_conversion)

            btn_layout.addWidget(self.run_btn)
            main_layout.addLayout(btn_layout)

        def browse_file(self):
            file_dialog = QFileDialog()
            file_dialog.setNameFilter("EBSD Files (*.ctf *.ang);;All Files (*)")
            if file_dialog.exec():
                files = file_dialog.selectedFiles()
                if files:
                    self.file_edit.setText(files[0])

        def log(self, message):
            self.log_text.append(message)
            self.log_text.verticalScrollBar().setValue(self.log_text.verticalScrollBar().maximum())

        def run_conversion(self):
            filename = self.file_edit.text().strip()
            if not filename:
                QMessageBox.warning(self, "Missing File", "Please select an input EBSD file.")
                return

            # Create args object
            class Args:
                pass

            args = Args()
            args.filename = filename
            args.frequency = self.freq_spin.value()
            args.thickness = self.thick_spin.value()
            args.nx_min = self.nx_min_spin.value()
            args.nx_max = self.nx_max_spin.value()
            args.ny_min = self.ny_min_spin.value()
            args.ny_max = self.ny_max_spin.value()
            args.N_grain_sample = self.n_grain_spin.value()
            args.UMAT = self.umat_check.isChecked()
            args.aster = self.aster_check.isChecked()
            args.multiphase = self.multiphase_check.isChecked()
            args.stereographic_projection = self.stereo_check.isChecked()
            args.odf = self.odf_check.isChecked()
            args.neper = self.neper_check.isChecked()

            # Disable run button
            self.run_btn.setEnabled(False)

            # Run in worker thread
            self.worker = ConversionWorker(args)
            self.worker.progress.connect(self.log)
            self.worker.finished.connect(self.conversion_finished)
            self.worker.start()
            
        def conversion_finished(self):
            self.run_btn.setEnabled(True)
            
    app = QApplication(sys.argv)
    window = EBSD2MOOSEGUI()
    window.show()
    sys.exit(app.exec())

parser = argparse.ArgumentParser(prog='EBSD2MOOSE',
                                 description='Convert EBSD ctf and ang files to MOOSE or UMAT Euler angles files',
                                 epilog='Nicolò Grilli, University of Bristol')

parser.add_argument('filename', nargs='?', help='Input EBSD file (ctf or ang format)')
parser.add_argument('-f', '--frequency', type=int, default=1, help='Frequency parameter')
parser.add_argument('-t', '--thickness', type=int, default=1, help='Thickness parameter')
parser.add_argument('-nx_min', '--nx_min', type=int, default=-1, help='Minimum nx')
parser.add_argument('-nx_max', '--nx_max', type=int, default=-1, help='Maximum nx')
parser.add_argument('-ny_min', '--ny_min', type=int, default=-1, help='Minimum ny')
parser.add_argument('-ny_max', '--ny_max', type=int, default=-1, help='Maximum ny')
parser.add_argument('-N_grain_sample', '--N_grain_sample', type=int, default=1, help='Number of grains for ODF sampling')
parser.add_argument('-UMAT', '--UMAT', action='store_true', help='Generate UMAT format output')
parser.add_argument('-aster', '--aster', action='store_true', help='Generate Aster format')
parser.add_argument('-multiphase', '--multiphase', action='store_true', help='Generate multiphase file')
parser.add_argument('-stereographic_projection', '--stereographic_projection', action='store_true', help='Generate stereographic projection')
parser.add_argument('-odf', '--odf', action='store_true', help='Generate ODF')
parser.add_argument('-neper', '--neper', action='store_true', help='Use Neper mesh file')
parser.add_argument('--gui', action='store_true', help='Launch graphical interface')

args = parser.parse_args()

# If --gui flag is provided or no filename is given, launch GUI
if args.gui or args.filename is None:
    main_gui()
else:
    run_conversion(args)
