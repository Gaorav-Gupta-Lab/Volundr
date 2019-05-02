#!/usr/bin/env python3
"""
Synthetic_Lethal_App.py v0.1.0
    Aug. 5, 2018
    Dennis A. Simpson
    This is a functional beta version.  Will work as expected.  Still needs cleaner GUI.
    The Target Search panel output still needs validation.

Code for a GUI to run Völundr Synthetic Lethal.
@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import datetime
import os
import collections
from contextlib import suppress
import dill
import wx
import wx.adv
from wx.lib.wordwrap import wordwrap
import wx.lib.sized_controls as sized_controls
import wx.lib.masked as masked
import wx.lib.intctrl
import subprocess

__author__ = 'Dennis A. Simpson'
__version__ = '0.3.0'
__package__ = 'Völundr'
__copyright__ = '(C) 2019'


class IntBoxes:
    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_int_caller(self, name, value=None):
        int_ctrl = wx.lib.intctrl.IntCtrl(self.main_panel, wx.ID_ANY, name=name, value=value, allow_none=1, min=0)
        int_ctrl.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        int_ctrl.SetColors(default_color=wx.BLACK, oob_color=wx.RED)
        return int_ctrl


class DataBoxes:
    """
    Generic ComboBox generator for project.  Ensures all boxes will look the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_box(self, data_list, name, value=None):
        if value is None:
            value = ""
        my_combo = wx.ComboBox(self.main_panel, wx.ID_DEFAULT, value, choices=data_list, style=wx.CB_SORT, name=name)
        my_combo.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        return my_combo


class MaskedDataBoxes:
    """
    Generic Masked ComboBox generator for project.  Ensures all boxes will look the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_masked_box(self, data_list, name, value=None):
        my_masked_combo = \
            masked.ComboBox(self.main_panel, wx.ID_DEFAULT, choices=data_list, style=wx.CB_SORT, name=name)
        my_masked_combo.SetCtrlParameters(invalidBackgroundColour="red", defaultValue=value, choiceRequired=True)
        my_masked_combo.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        return my_masked_combo


class ButtonGenerator:
    """
    Generic Button Generator for Project Insuring all Buttons are the same.
    """

    def __init__(self, main_panel):
        self.main_panel = main_panel

    def my_buttons(self, label):
        my_button = wx.Button(self.main_panel, wx.ID_ANY, label=label)
        my_button.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        return my_button


class BooleanValidator(wx.Validator):
    def __init__(self, parent, name):
        super(BooleanValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        return BooleanValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if "/" in text:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} must be a full path statement".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class PathValidator(wx.Validator):
    def __init__(self, parent, name):
        super(PathValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return PathValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if "/" in text:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} requires full path statements".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class TextValidator(wx.Validator):
    def __init__(self, parent, name):
        super(TextValidator, self).__init__()
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return TextValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if len(text) > 0:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} cannot be Null".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
        return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class NumValidator(wx.lib.intctrl.IntValidator):
    def __init__(self, parent, name):
        wx.lib.intctrl.IntValidator.__init__(self)
        self.parent = parent
        self.name = name

    def Clone(self):
        """ """
        return NumValidator(self.parent, self.name)

    def Validate(self, window):
        """ """
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()

        if text >= 0:
            textCtrl.SetBackgroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        else:
            message = "{} takes positive integers only".format(self.name)
            caption = "Invalid Input"
            dlg = wx.GenericMessageDialog(self.parent, message, caption, style=wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Show()
            dlg.Destroy()

            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()

            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class CommonControls:
    def __init__(self, parent):
        self.parent = parent
        self.default_data = wx.GetTopLevelParent(parent).default_data
        self.dataframe = wx.GetTopLevelParent(parent).dataframe
        self.button_generator = ButtonGenerator(parent)
        self.databoxes = DataBoxes(parent)
        self.my_masked_databoxes = MaskedDataBoxes(parent)
        self.my_int_box = IntBoxes(parent)

    def folder_selector(self, name, tip=None):
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name)
        input_ctrl.SetValidator(PathValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        working_dir_btn = self.button_generator.my_buttons(name)
        working_dir_btn.Bind(wx.EVT_BUTTON, lambda event, temp=working_dir_btn: self.select_folder(event, input_ctrl))

        return working_dir_btn, input_ctrl

    def file_selector(self, name, tip=None):
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name)
        input_ctrl.SetValidator(PathValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        fastq_file_btn = self.button_generator.my_buttons(name)
        fastq_file_btn.Bind(wx.EVT_BUTTON, lambda event, temp=fastq_file_btn: self.select_file(event, input_ctrl))

        return fastq_file_btn, input_ctrl

    def restricted_selector(self, name, value=None, tip=None):
        """
        Define the verbosity level for the logger.
        :return:
        """
        input_ctrl = self.my_masked_databoxes.my_masked_box(self.default_data[name], name, value)
        if tip is not None:
            input_ctrl.SetToolTip(tip)
        label = wx.StaticText(self.parent, wx.ID_ANY, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def add_line(self):
        label = wx.StaticText(self.parent, wx.ID_DEFAULT, "", style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL))

        label.SetMinSize((200, 10))
        line = wx.StaticLine(self.parent, wx.ID_ANY, style=wx.LI_HORIZONTAL, name="line")
        line.SetBackgroundColour('#E5E8E8')

        return label, line

    def text_control(self, name, tip=None):
        try:
            value = self.default_data[name]
        except KeyError:
            value = ""
        input_ctrl = self.databoxes.my_box(self.dataframe[name], name, value)
        input_ctrl.SetValidator(TextValidator(self.parent, input_ctrl.GetName()))
        if tip is not None:
            input_ctrl.SetToolTip(tip)

        label = wx.StaticText(self.parent, wx.ID_DEFAULT, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))
        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def int_control(self, name, tip=None):
        input_ctrl = self.my_int_box.my_int_caller(name, self.default_data[name])
        if tip is not None:
            input_ctrl.SetToolTip(tip)

        label = wx.StaticText(self.parent, wx.ID_DEFAULT, name.strip("--"), style=wx.ALIGN_RIGHT)
        label.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.NORMAL, faceName='Inconsolata'))

        label.SetForegroundColour("blue")
        label.SetMinSize((200, 20))

        return label, input_ctrl

    def select_folder(self, event, ctrl):
        dlg = wx.DirDialog(self.parent, "Choose a Folder")
        if dlg.ShowModal() == wx.ID_OK:
            ctrl.SetValue(dlg.GetPath())
            dlg.Destroy()

    def select_file(self, event, ctrl):
        dlg = wx.FileDialog(self.parent, "Choose a file")
        if dlg.ShowModal() == wx.ID_OK:
            ctrl.SetValue(dlg.GetPath())
            dlg.Destroy()

    def int_control_builder(self, name, default_value):
        name = name
        try:
            value = default_value[name]
        except KeyError:
            value = None

        input_ctrl = self.my_int_box.my_int_caller(name, value)

        return input_ctrl


class StatisticsPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(StatisticsPanel, self).__init__(parent, wx.ID_ANY, name="StatisticsPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Statistics Parameters")

        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list =\
            [my_controls.file_selector("Völundr_API", tip="Full Path Volundr.py"),
             my_controls.folder_selector("--Options_File",
                                         tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.file_selector("--Target_File", tip="Full Path to Target File"),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.file_selector("--Master_Index_File", tip="Full Path to Master Index File"),
             my_controls.file_selector("--Index_File", tip="Full Path to Multiplex Index.bed File"),
             my_controls.add_line(),
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name", tip="Name for run. For Statistics this must match Job Name on Count"
                                                        " Files"),
             my_controls.int_control("--Spawn", tip="How many CPU's or Threads.  Leave at 1 for Statistics"),
             my_controls.restricted_selector("--Species", tip="Currently must be Human or Mouse"),
             my_controls.text_control("--Control_Sample", tip="Name of biological control in Index File"),
             my_controls.int_control("--Target_Mismatch", tip="Number of mismatches allowed in target sequence"),
             my_controls.text_control("--Target_Length",
                                      tip="Length of target from Target_File.  Accepted values are an integer or "
                                          "'Variable'.  Default '20'"),
             my_controls.int_control("--Target_Start",
                                     tip="Location of target sequence in Target File sequences. Not used if "
                                         "Target_Length=Variable.  Default=20"),
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))
            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)
            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class TargetSearchPanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(TargetSearchPanel, self).__init__(parent, wx.ID_ANY, name="TargetSearchPanel")
        panel_sizer = wx.BoxSizer(wx.VERTICAL)
        title_sizer = wx.BoxSizer(wx.HORIZONTAL)
        panel_grid_sizer = wx.GridBagSizer(0, 0)
        my_controls = CommonControls(self)

        title_bmp = wx.ArtProvider.GetBitmap(wx.ART_INFORMATION, wx.ART_OTHER, (32, 32))
        title_icon = wx.StaticBitmap(self, wx.ID_ANY, title_bmp)
        title = wx.StaticText(self, wx.ID_ANY, "Target Search Parameters")
        title.SetForegroundColour("blue")
        title.SetSize(title.GetBestSize())
        title_sizer.Add(title_icon, 0, wx.ALL)
        title_sizer.Add(title, 0, wx.EXPAND)

        panel_sizer.Add(title_sizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        panel_sizer.Add(wx.StaticLine(self, ), 0, wx.ALL | wx.EXPAND)

        # Build a list of our control objects that create the widgets.  The order on the form is the same as this order
        widget_build_list =\
            [my_controls.file_selector("Völundr_API", tip="Full Path to Volundr.py"),
             my_controls.folder_selector("--Options_File", tip="Location of this file.  No trailing slash. Generally full path is needed."),
             my_controls.folder_selector("--Working_Folder", tip="Full Path to Working Folder, no Trailing Slash"),
             my_controls.file_selector("--FASTQ1", tip="Full Path multiplexed FASTQ file"),
             my_controls.file_selector("--Target_File", tip="Full Path to Target File"),
             my_controls.file_selector("--Master_Index_File", tip="Full Path to Master Index File"),
             my_controls.file_selector("--Index_File", tip="Full Path to Multiplex Index.bed File"),
             my_controls.add_line(),
             my_controls.restricted_selector("--Verbose", "INFO"),
             my_controls.text_control("--Job_Name", tip="Name for run.  No special characters"),
             my_controls.int_control("--Spawn", tip="How many parallel jobs.  Maximum is n-1 CPU's or Threads."),
             my_controls.restricted_selector("--Species", tip="Currently must be Human or Mouse"),
             my_controls.restricted_selector("--Analyze_Unknowns", tip="Process reads with unidentifiable indices?"),
             my_controls.restricted_selector("--Delete_Demultiplexed_FASTQ", "True", tip="Delete the temporary demultiplexed FASTQ Files.  Default is True"),
             my_controls.restricted_selector("--Compress", tip="If keeping demultiplexe FASTQ files should they be compressed?"),
             my_controls.restricted_selector("--RevComp", tip="Are the sequence reads the reverse compliment of the sgGuide Targets?"),
             my_controls.add_line(),
             my_controls.int_control("--Target_Mismatch", tip="Number of mismatches allowed in target sequence"),
             my_controls.int_control("--Min_Length", tip="Minimum length of sequence read."),
             my_controls.text_control("--Target_Length", tip="Length of target from Target_File.  Accepted values are an integer or 'Variable'.  Default '20'"),
             my_controls.int_control("--Target_Start", tip="Location of target sequence in Target File sequences. Not used if Target_Length=Variable.  Default=20"),
             my_controls.int_control("--Index_Mismatch", tip="Number of mismatches allowed in library indices.  Default=1"),
             my_controls.int_control("--Target_Padding", tip="Additional nucleotides to include in unknown sequence if AnchorSeq is not found."),
             my_controls.int_control("--Expected_Position", tip="Location in the sequence read the target is expected.  Used if AnchorSeq is not found."),
             my_controls.text_control("--AnchorSeq", tip="Nucleotides 5' of target used in search.  Default 'AAACACCG'"),
             my_controls.int_control("--AnchorMismatch", tip="Number of mismatches allowed in anchor sequence.  Default=1"),
             my_controls.int_control("--AnchorStart", tip="Position in sequence read to being search for Anchor Sequence"),
             my_controls.int_control("--AnchorStop", tip="Position in sequence read to end search for Anchor Sequence"),
             ]

        # Put the widgets on the form
        self.control_dict = collections.defaultdict(tuple)
        static_line_sizer = wx.BoxSizer(wx.HORIZONTAL)
        static_line_sizer.Add(wx.StaticLine(self), 0, wx.EXPAND)
        for i in range(len(widget_build_list)):
            widget0 = widget_build_list[i][0]
            widget1 = widget_build_list[i][1]
            panel_grid_sizer.Add(width=20, height=0, pos=(i, 0))

            panel_grid_sizer.Add(widget0, pos=(i, 1), flag=wx.ALL | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=3)
            panel_grid_sizer.Add(widget1, pos=(i, 2),  flag=wx.EXPAND | wx.ALL, border=3)

            panel_grid_sizer.Add(width=10, height=0, pos=(i, 3))
            self.control_dict[i] = (widget0, widget1)

        # panel_grid_sizer.AddGrowableCol(1)
        panel_grid_sizer.AddGrowableCol(2)
        panel_sizer.Add(panel_grid_sizer, 0, wx.ALL | wx.EXPAND)
        self.SetSizerAndFit(panel_sizer)


class WelcomePanel(sized_controls.SizedScrolledPanel):
    def __init__(self, parent):
        super(WelcomePanel, self).__init__(parent, wx.ID_ANY, name="WelcomePanel")
        self.SetBackgroundStyle(wx.BG_STYLE_ERASE)
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.add_background)
        self.SetSizerAndFit(self.box)

    def add_background(self, event):
        """
         Add a picture to the background
        """
        dc = wx.ClientDC(self)
        rect = self.GetUpdateRegion().GetBox()
        dc.SetClippingRegion(rect)

        os.chdir('..')
        # bkgrnd_image = "{}/img/Volundr_Splash.png".format(os.getcwd())
        bkgrnd_image = "/OneDrive_UNC/Projects/Volundr/img/Volundr_Splash.png"

        # If the background file is missing the installation is likely damaged.
        if not os.path.isfile(bkgrnd_image):
            print("Background Image File {} not Found.  Installation Not Correct".format(bkgrnd_image))
            raise SystemExit(1)

        bmp = wx.Bitmap(bkgrnd_image)
        cliWidth, cliHeight = self.GetClientSize()

        bmpWide = bmp.GetWidth()
        bmpHeight = bmp.GetHeight()
        img = bmp.ConvertToImage()
        bmpWide = int((cliWidth / bmpWide) * bmpWide)
        bmpHeight = int((cliHeight / bmpHeight) * bmpHeight)

        bmp = wx.Bitmap(img.Scale(bmpWide, bmpHeight))

        xPos = (cliWidth - (bmpWide)) / 2
        yPos = (cliHeight - (bmpHeight)) / 2
        dc.DrawBitmap(bmp, xPos, yPos)


class MainFrame(wx.Frame):
    def __init__(self):
        self.dir_path = os.path.dirname(os.path.abspath(__file__))
        self.pickle_file = "{0}{1}pickles{1}CRISPR_parameters.pkl".format(self.dir_path, os.sep)

        # Setup a frame that is centered on the screen and opens at 75% of full screen.
        window_name = "Volundr Synthetic Lethal GUI  v{}".format(__version__)
        display_x, display_y = wx.DisplaySize()
        scale = 0.75

        wx.Frame.__init__(self, None, wx.ID_ANY, title=window_name, size=(display_x*scale, display_y*scale))
        self.Centre()
        self.build_menu_bar()

        self.dataframe = self.dataframe_build()
        self.default_data = self.default_dict_build()
        print(self.GetSize(), wx.DisplaySize(), wx.version())

        self.welcome_panel = WelcomePanel(self)
        self.target_search_panel = TargetSearchPanel(self)
        self.statistics_panel = StatisticsPanel(self)

        self.target_search_panel.Hide()
        self.statistics_panel.Hide()

        self.panel_dict = \
            {"TargetSearchPanel": self.target_search_panel,
             "WelcomePanel": self.welcome_panel,
             "StatisticsPanel": self.statistics_panel}

        self.main_sizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = self.welcome_panel
        self.main_sizer.Add(self.panel, 1, wx.EXPAND)

        self.SetSizer(self.main_sizer)

    @staticmethod
    def default_dict_build():
        """
        This builds a dictionary with default values that are library type and version specific.
        :return:
        """
        version1_default_dict = \
            {"--AnchorStart": 100, "--AnchorStop": 125, "--Expected_Position": 112, "--Target_Start": 20,
             "--Target_Length": "20", "--AnchorSeq": "AAACACCG", "--Delete_Demultiplexed_FASTQ": ["True", "False"],
             "--Species": ["Human", "Mouse"], "--Analyze_Unknowns": ["True", "False"], "--Target_Mismatch": 1,
             "--Verbose": ["INFO", "DEBUG", "ERROR", "WARN"], "--Spawn": 1, "--Compress": ["True", "False"],
             "--Index_Mismatch": 1, "--Target_Padding": 2, "--AnchorMismatch": 1, "--RevComp": ["True", "False"],
             "--Min_Length": 120
             }
        return version1_default_dict

    def dataframe_build(self):
        try:
            with open(self.pickle_file, 'rb') as file:
                dataframe_dict = dill.load(file)

        except FileNotFoundError:
            column_names = \
                ["--FASTQ1", "--Index_File", "--Target_File", "--Master_Index_File", "--Working_Folder", "--Verbose",
                 "--Job_Name", "--Spawn", "--Target_Search", "--Analyze_Unknowns", "--Statistics", "--RevComp"
                 "--Combine_Replicates", "--Species", "--Control_Sample", "--Min_Length", "Völundr_API",
                 "--Target_Length", "--Target_Start", "--Index_Mismatch", "--Target_Padding", "--Target_Mismatch",
                 "--Expected_Position", "--AnchorSeq", "--AnchorMismatch", "--AnchorStart", "--AnchorStop"]

            dataframe_dict = collections.defaultdict(list)
            for key in column_names:
                dataframe_dict[key] = []
        return dataframe_dict

    def build_menu_bar(self):
        menu_bar = wx.MenuBar()

        # Build the "File" menu.
        file_menu = wx.Menu()
        welcome_app = file_menu.Append(wx.ID_ANY, "Home", "Go to Welcome Screen")
        save_app = file_menu.Append(wx.ID_SAVE, "&Save\tCtrl-S", "Display a hint about this function")
        run_app = file_menu.Append(wx.ID_EXECUTE, "&Run\tCtrl-R", "Display a hint about this function")
        file_menu.AppendSeparator()
        exit_app = file_menu.Append(wx.ID_EXIT, "E&xit\tAlt-x", "Close window and exit program.")
        menu_bar.Append(file_menu, "&File")

        # Build the "Tools" Menu
        tools_menu = wx.Menu()
        target_search_app = tools_menu.Append(wx.ID_ANY, "TargetSearch", "Search for CRISPR Guides")
        analyze_counts_app = tools_menu.Append(wx.ID_ANY, "AnalyzeCounts", "Do Some Statistics on Counts")
        menu_bar.Append(tools_menu, "&Tools")

        # Build the "Help" menu
        help_menu = wx.Menu()
        about_item = help_menu.Append(wx.ID_ABOUT)
        menu_bar.Append(help_menu, "&Help")
        self.SetMenuBar(menu_bar)

        # Bind Menu Items to Actions
        self.Bind(wx.EVT_MENU, self.save_parameter_file, save_app)
        self.Bind(wx.EVT_MENU, self.on_run, run_app)
        self.Bind(wx.EVT_MENU, self.on_exit, exit_app)
        self.Bind(wx.EVT_MENU, self.on_about, about_item)

        # Panel switching
        self.Bind(wx.EVT_MENU, lambda event, temp=target_search_app:
                  self.switch_panels(event, self.target_search_panel.GetName()), target_search_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=analyze_counts_app:
                  self.switch_panels(event, self.statistics_panel.GetName()), analyze_counts_app)
        self.Bind(wx.EVT_MENU, lambda event, temp=welcome_app:
                  self.switch_panels(event, self.welcome_panel.GetName()), welcome_app)

    def switch_panels(self, event, switch_id):
        self.main_sizer.Detach(self.panel)
        self.panel.Hide()
        self.panel = self.panel_dict[switch_id]
        self.main_sizer.Add(self.panel, 1, wx.EXPAND)
        self.panel.Show()
        self.panel.Fit()
        self.Layout()

    def on_exit(self, event):
        """Close the frame, terminating the application."""
        self.Close(True)

    def on_about(self, event):
        about_info = wx.adv.AboutDialogInfo()

        about_info.SetName("Synthetic Lethal App")
        about_info.SetCopyright(__copyright__)

        about_info.SetDescription(wordwrap(
            "GUI for Synthetic Lethal Module of Volundr.\nVersion: {}".format(__version__), 350, wx.ClientDC(self)))
        about_info.SetDevelopers([__author__])
        # about_info.License = wordwrap("Completely and totally open source!", 500, wx.ClientDC(self))
        wx.adv.AboutBox(about_info)

    def on_run(self, event):
        shellfile_name = self.save_parameter_file("Run")
        print(shellfile_name)
        if shellfile_name is not None:
            subprocess.run([shellfile_name], shell=True)

    def save_parameter_file(self, event):
        def outstring_build():
            file_body = "{}{}\n".format(target_search, analyze_counts)
            working_folder = ""
            api = ""
            job = ""
            options = ""

            for i in range(len(panel.control_dict)):
                d = panel.control_dict[i][1]

                # Validate data.  If validation fails this will take us back to the form.
                with suppress(AttributeError):
                    if not d.GetValidator().Validate(d):
                        return job, working_folder, file_body, False

                ctrl_name = d.GetName()
                if ctrl_name == 'line':
                    continue
                ctrl_value = str(d.GetValue()).strip()

                # Linux does not add the trailing slash to folders.  I prefer to have one.
                if ctrl_name == "Völundr_API":
                    api = ctrl_value
                elif ctrl_name == "--Working_Folder":
                    working_folder = ctrl_value
                    ctrl_value += "/"
                elif ctrl_name == "--Job_Name":
                    job = ctrl_value
                elif ctrl_name == "--Options_File":
                    options = ctrl_value

                if not ctrl_name == "Völundr_API" and not ctrl_name == "--Options_File":
                    file_body += "{}\t{}\n".format(ctrl_name, ctrl_value)
                self.dataframe[ctrl_name].append(d.GetValue())
                self.dataframe[ctrl_name] = list(set(self.dataframe[ctrl_name]))

            return job, working_folder, file_body, True, api, options

        if self.target_search_panel.IsShown():
            submodule = "Target_Search"
            panel = self.target_search_panel
            target_search = "--Target_Search\tTrue\n"
            analyze_counts = "--Statistics\tFalse\n"

        elif self.statistics_panel.IsShown():
            submodule = "Analyze_Counts"
            panel = self.statistics_panel
            target_search = "--Target_Search\tFalse\n"
            analyze_counts = "--Statistics\tTrue\n"

        job_name, working_dir, parameter_body, validation_pass, api_path, options_path = outstring_build()

        # If the validation fails take the user back to the panel so they can correct the error(s)
        if not validation_pass:
            return

        date = datetime.datetime.now().strftime("%Y%m%d")
        outfile_name = "{}{}run_{}_{}.sh".format(options_path, os.sep, job_name, date)
        shebang = "#!/bin/bash\n" \
                  "#Parameter file to run Völundr Synthetic Lethal module {}\n" \
                  "#File generated {}\n\n".format(submodule, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        cmd = "python3 {} --options_file {}\nexit\n\n".format(api_path, outfile_name)

        # Save the current choices in the pickles folder.
        with open(self.pickle_file, 'wb') as file:
            dill.dump(self.dataframe, file, protocol=-1)

        # Write the parameter file.
        outstring = "{}{}{}".format(shebang, cmd, parameter_body)
        filename = "{}_{}.sh".format(job_name, date)

        if event == "Run" and validation_pass:
            # parameter_file = open("{}{}{}".format(working_dir, os.sep, filename), 'w')
            parameter_file = open(outfile_name, 'w')
            parameter_file.write(outstring)
            parameter_file.close()
            return outfile_name

        with wx.FileDialog(self, "Save Völundr Options File",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT, defaultFile=outfile_name) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return  # the user changed their mind

            # save the current contents in the file
            pathname = fileDialog.GetPath()

            try:
                parameter_file = open("{}".format(pathname), 'w')
                parameter_file.write(outstring)
                parameter_file.close()
            except IOError:
                wx.LogError("Cannot save {} in {}".format(filename, pathname))


def main():
    app = wx.App()
    frame = MainFrame()
    frame.Show()
    app.MainLoop()
    # wx.lib.inspection.InspectionTool().Show()


if __name__ == '__main__':
    main()
