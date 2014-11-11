namespace DTA_Mosaicing_GUI
{
    partial class SettingsGUI
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.buttonGo = new System.Windows.Forms.Button();
            this.buttonExit = new System.Windows.Forms.Button();
            this.tabControlMain = new System.Windows.Forms.TabControl();
            this.tabVidSource = new System.Windows.Forms.TabPage();
            this.groupBoxOutput = new System.Windows.Forms.GroupBox();
            this.buttonBrowseSaveIm = new System.Windows.Forms.Button();
            this.textBoxSaveMosaic = new System.Windows.Forms.TextBox();
            this.radioButtonSaveMosaic = new System.Windows.Forms.RadioButton();
            this.radioButtonDisplayOnly = new System.Windows.Forms.RadioButton();
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.buttonBrowseDir = new System.Windows.Forms.Button();
            this.buttonBrowseVidFile = new System.Windows.Forms.Button();
            this.textBoxDir = new System.Windows.Forms.TextBox();
            this.textBoxVidFile = new System.Windows.Forms.TextBox();
            this.radioButtonDir = new System.Windows.Forms.RadioButton();
            this.radioButtonVideo = new System.Windows.Forms.RadioButton();
            this.radioButtonCam = new System.Windows.Forms.RadioButton();
            this.tabTransEstimator = new System.Windows.Forms.TabPage();
            this.groupBox4 = new System.Windows.Forms.GroupBox();
            this.label4 = new System.Windows.Forms.Label();
            this.textBoxCorrQuality = new System.Windows.Forms.TextBox();
            this.groupBox3 = new System.Windows.Forms.GroupBox();
            this.labelBlobScales = new System.Windows.Forms.Label();
            this.textBoxBlobScales = new System.Windows.Forms.TextBox();
            this.label6 = new System.Windows.Forms.Label();
            this.comboBoxPatchSize = new System.Windows.Forms.ComboBox();
            this.label5 = new System.Windows.Forms.Label();
            this.comboBoxDescriptor = new System.Windows.Forms.ComboBox();
            this.label3 = new System.Windows.Forms.Label();
            this.textBoxNumCorners = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.comboBoxFeatureDetector = new System.Windows.Forms.ComboBox();
            this.groupBox2 = new System.Windows.Forms.GroupBox();
            this.labelLMiters = new System.Windows.Forms.Label();
            this.textBoxLMIterations = new System.Windows.Forms.TextBox();
            this.checkBoxLM = new System.Windows.Forms.CheckBox();
            this.comboBoxTransType = new System.Windows.Forms.ComboBox();
            this.label1 = new System.Windows.Forms.Label();
            this.tabRenderer = new System.Windows.Forms.TabPage();
            this.groupBox5 = new System.Windows.Forms.GroupBox();
            this.label13 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.textBoxFullFrameUpdateFreq = new System.Windows.Forms.TextBox();
            this.checkBoxSkipFrames = new System.Windows.Forms.CheckBox();
            this.label7 = new System.Windows.Forms.Label();
            this.comboBoxWarpMethod = new System.Windows.Forms.ComboBox();
            this.labelFeatherRad = new System.Windows.Forms.Label();
            this.textBoxFeatherRad = new System.Windows.Forms.TextBox();
            this.labelDijkstraScale = new System.Windows.Forms.Label();
            this.textBoxDijkstraScale = new System.Windows.Forms.TextBox();
            this.label11 = new System.Windows.Forms.Label();
            this.textBoxMosaicY = new System.Windows.Forms.TextBox();
            this.label10 = new System.Windows.Forms.Label();
            this.textBoxMosaicX = new System.Windows.Forms.TextBox();
            this.label9 = new System.Windows.Forms.Label();
            this.textBoxMaxFrames = new System.Windows.Forms.TextBox();
            this.checkBoxIncrementalRenderer = new System.Windows.Forms.CheckBox();
            this.label8 = new System.Windows.Forms.Label();
            this.comboBoxRenderer = new System.Windows.Forms.ComboBox();
            this.tabPageEval = new System.Windows.Forms.TabPage();
            this.groupBox6 = new System.Windows.Forms.GroupBox();
            this.label18 = new System.Windows.Forms.Label();
            this.comboBoxEvalFunction = new System.Windows.Forms.ComboBox();
            this.groupBox7 = new System.Windows.Forms.GroupBox();
            this.textBoxMaxRansacIters = new System.Windows.Forms.TextBox();
            this.textBoxMaxSearchForTransform = new System.Windows.Forms.TextBox();
            this.label14 = new System.Windows.Forms.Label();
            this.label15 = new System.Windows.Forms.Label();
            this.tabControlMain.SuspendLayout();
            this.tabVidSource.SuspendLayout();
            this.groupBoxOutput.SuspendLayout();
            this.groupBox1.SuspendLayout();
            this.tabTransEstimator.SuspendLayout();
            this.groupBox4.SuspendLayout();
            this.groupBox3.SuspendLayout();
            this.groupBox2.SuspendLayout();
            this.tabRenderer.SuspendLayout();
            this.groupBox5.SuspendLayout();
            this.tabPageEval.SuspendLayout();
            this.groupBox6.SuspendLayout();
            this.groupBox7.SuspendLayout();
            this.SuspendLayout();
            // 
            // buttonGo
            // 
            this.buttonGo.DialogResult = System.Windows.Forms.DialogResult.OK;
            this.buttonGo.Location = new System.Drawing.Point(424, 592);
            this.buttonGo.Name = "buttonGo";
            this.buttonGo.Size = new System.Drawing.Size(104, 32);
            this.buttonGo.TabIndex = 0;
            this.buttonGo.Text = "&Start capture";
            this.buttonGo.UseVisualStyleBackColor = true;
            this.buttonGo.Click += new System.EventHandler(this.buttonGo_Click);
            // 
            // buttonExit
            // 
            this.buttonExit.DialogResult = System.Windows.Forms.DialogResult.Cancel;
            this.buttonExit.Location = new System.Drawing.Point(544, 592);
            this.buttonExit.Name = "buttonExit";
            this.buttonExit.Size = new System.Drawing.Size(96, 32);
            this.buttonExit.TabIndex = 1;
            this.buttonExit.Text = "E&xit";
            this.buttonExit.UseVisualStyleBackColor = true;
            this.buttonExit.Click += new System.EventHandler(this.buttonExit_Click);
            // 
            // tabControlMain
            // 
            this.tabControlMain.Controls.Add(this.tabVidSource);
            this.tabControlMain.Controls.Add(this.tabTransEstimator);
            this.tabControlMain.Controls.Add(this.tabRenderer);
            this.tabControlMain.Controls.Add(this.tabPageEval);
            this.tabControlMain.Location = new System.Drawing.Point(8, 8);
            this.tabControlMain.Name = "tabControlMain";
            this.tabControlMain.SelectedIndex = 0;
            this.tabControlMain.Size = new System.Drawing.Size(632, 576);
            this.tabControlMain.TabIndex = 2;
            // 
            // tabVidSource
            // 
            this.tabVidSource.Controls.Add(this.groupBoxOutput);
            this.tabVidSource.Controls.Add(this.groupBox1);
            this.tabVidSource.Location = new System.Drawing.Point(4, 22);
            this.tabVidSource.Name = "tabVidSource";
            this.tabVidSource.Padding = new System.Windows.Forms.Padding(3);
            this.tabVidSource.Size = new System.Drawing.Size(624, 550);
            this.tabVidSource.TabIndex = 0;
            this.tabVidSource.Text = "Video source";
            this.tabVidSource.UseVisualStyleBackColor = true;
            // 
            // groupBoxOutput
            // 
            this.groupBoxOutput.Controls.Add(this.buttonBrowseSaveIm);
            this.groupBoxOutput.Controls.Add(this.textBoxSaveMosaic);
            this.groupBoxOutput.Controls.Add(this.radioButtonSaveMosaic);
            this.groupBoxOutput.Controls.Add(this.radioButtonDisplayOnly);
            this.groupBoxOutput.Location = new System.Drawing.Point(8, 176);
            this.groupBoxOutput.Name = "groupBoxOutput";
            this.groupBoxOutput.Size = new System.Drawing.Size(608, 120);
            this.groupBoxOutput.TabIndex = 10;
            this.groupBoxOutput.TabStop = false;
            this.groupBoxOutput.Text = "Mosaic output";
            // 
            // buttonBrowseSaveIm
            // 
            this.buttonBrowseSaveIm.Location = new System.Drawing.Point(496, 64);
            this.buttonBrowseSaveIm.Name = "buttonBrowseSaveIm";
            this.buttonBrowseSaveIm.Size = new System.Drawing.Size(72, 20);
            this.buttonBrowseSaveIm.TabIndex = 9;
            this.buttonBrowseSaveIm.Text = "Browse...";
            this.buttonBrowseSaveIm.UseVisualStyleBackColor = true;
            this.buttonBrowseSaveIm.Click += new System.EventHandler(this.buttonBrowseSaveIm_Click);
            // 
            // textBoxSaveMosaic
            // 
            this.textBoxSaveMosaic.Location = new System.Drawing.Point(152, 64);
            this.textBoxSaveMosaic.Name = "textBoxSaveMosaic";
            this.textBoxSaveMosaic.Size = new System.Drawing.Size(344, 20);
            this.textBoxSaveMosaic.TabIndex = 7;
            this.textBoxSaveMosaic.TextChanged += new System.EventHandler(this.textBoxSaveMosaic_TextChanged);
            // 
            // radioButtonSaveMosaic
            // 
            this.radioButtonSaveMosaic.AutoSize = true;
            this.radioButtonSaveMosaic.Location = new System.Drawing.Point(24, 64);
            this.radioButtonSaveMosaic.Name = "radioButtonSaveMosaic";
            this.radioButtonSaveMosaic.Size = new System.Drawing.Size(124, 17);
            this.radioButtonSaveMosaic.TabIndex = 5;
            this.radioButtonSaveMosaic.TabStop = true;
            this.radioButtonSaveMosaic.Text = "&Save to this directory";
            this.radioButtonSaveMosaic.UseVisualStyleBackColor = true;
            this.radioButtonSaveMosaic.CheckedChanged += new System.EventHandler(this.radioButtonSaveMosaic_CheckedChanged);
            // 
            // radioButtonDisplayOnly
            // 
            this.radioButtonDisplayOnly.AutoSize = true;
            this.radioButtonDisplayOnly.Checked = true;
            this.radioButtonDisplayOnly.Location = new System.Drawing.Point(24, 32);
            this.radioButtonDisplayOnly.Name = "radioButtonDisplayOnly";
            this.radioButtonDisplayOnly.Size = new System.Drawing.Size(126, 17);
            this.radioButtonDisplayOnly.TabIndex = 3;
            this.radioButtonDisplayOnly.TabStop = true;
            this.radioButtonDisplayOnly.Text = "Output to &display only";
            this.radioButtonDisplayOnly.UseVisualStyleBackColor = true;
            this.radioButtonDisplayOnly.CheckedChanged += new System.EventHandler(this.radioButtonDisplayOnly_CheckedChanged);
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.buttonBrowseDir);
            this.groupBox1.Controls.Add(this.buttonBrowseVidFile);
            this.groupBox1.Controls.Add(this.textBoxDir);
            this.groupBox1.Controls.Add(this.textBoxVidFile);
            this.groupBox1.Controls.Add(this.radioButtonDir);
            this.groupBox1.Controls.Add(this.radioButtonVideo);
            this.groupBox1.Controls.Add(this.radioButtonCam);
            this.groupBox1.Location = new System.Drawing.Point(8, 16);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Size = new System.Drawing.Size(608, 144);
            this.groupBox1.TabIndex = 3;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "Image source";
            // 
            // buttonBrowseDir
            // 
            this.buttonBrowseDir.Location = new System.Drawing.Point(496, 96);
            this.buttonBrowseDir.Name = "buttonBrowseDir";
            this.buttonBrowseDir.Size = new System.Drawing.Size(72, 20);
            this.buttonBrowseDir.TabIndex = 9;
            this.buttonBrowseDir.Text = "Browse...";
            this.buttonBrowseDir.UseVisualStyleBackColor = true;
            this.buttonBrowseDir.Click += new System.EventHandler(this.buttonBrowseDir_Click);
            // 
            // buttonBrowseVidFile
            // 
            this.buttonBrowseVidFile.Location = new System.Drawing.Point(496, 64);
            this.buttonBrowseVidFile.Name = "buttonBrowseVidFile";
            this.buttonBrowseVidFile.Size = new System.Drawing.Size(72, 20);
            this.buttonBrowseVidFile.TabIndex = 8;
            this.buttonBrowseVidFile.Text = "Browse...";
            this.buttonBrowseVidFile.UseVisualStyleBackColor = true;
            this.buttonBrowseVidFile.Click += new System.EventHandler(this.buttonBrowse_Click);
            // 
            // textBoxDir
            // 
            this.textBoxDir.Location = new System.Drawing.Point(152, 96);
            this.textBoxDir.Name = "textBoxDir";
            this.textBoxDir.Size = new System.Drawing.Size(344, 20);
            this.textBoxDir.TabIndex = 7;
            this.textBoxDir.TextChanged += new System.EventHandler(this.textBoxDir_TextChanged);
            // 
            // textBoxVidFile
            // 
            this.textBoxVidFile.Location = new System.Drawing.Point(152, 64);
            this.textBoxVidFile.Name = "textBoxVidFile";
            this.textBoxVidFile.Size = new System.Drawing.Size(344, 20);
            this.textBoxVidFile.TabIndex = 6;
            this.textBoxVidFile.TextChanged += new System.EventHandler(this.textBoxVidFile_TextChanged);
            // 
            // radioButtonDir
            // 
            this.radioButtonDir.AutoSize = true;
            this.radioButtonDir.Location = new System.Drawing.Point(24, 96);
            this.radioButtonDir.Name = "radioButtonDir";
            this.radioButtonDir.Size = new System.Drawing.Size(97, 17);
            this.radioButtonDir.TabIndex = 5;
            this.radioButtonDir.TabStop = true;
            this.radioButtonDir.Text = "&Image directory";
            this.radioButtonDir.UseVisualStyleBackColor = true;
            this.radioButtonDir.CheckedChanged += new System.EventHandler(this.radioButtonDir_CheckedChanged);
            // 
            // radioButtonVideo
            // 
            this.radioButtonVideo.AutoSize = true;
            this.radioButtonVideo.Location = new System.Drawing.Point(24, 64);
            this.radioButtonVideo.Name = "radioButtonVideo";
            this.radioButtonVideo.Size = new System.Drawing.Size(68, 17);
            this.radioButtonVideo.TabIndex = 4;
            this.radioButtonVideo.TabStop = true;
            this.radioButtonVideo.Text = "&Video file";
            this.radioButtonVideo.UseVisualStyleBackColor = true;
            this.radioButtonVideo.CheckedChanged += new System.EventHandler(this.radioButtonVideo_CheckedChanged);
            // 
            // radioButtonCam
            // 
            this.radioButtonCam.AutoSize = true;
            this.radioButtonCam.Checked = true;
            this.radioButtonCam.Location = new System.Drawing.Point(24, 32);
            this.radioButtonCam.Name = "radioButtonCam";
            this.radioButtonCam.Size = new System.Drawing.Size(106, 17);
            this.radioButtonCam.TabIndex = 3;
            this.radioButtonCam.TabStop = true;
            this.radioButtonCam.Text = "Attached &camera";
            this.radioButtonCam.UseVisualStyleBackColor = true;
            this.radioButtonCam.CheckedChanged += new System.EventHandler(this.radioButtonCam_CheckedChanged);
            // 
            // tabTransEstimator
            // 
            this.tabTransEstimator.Controls.Add(this.groupBox7);
            this.tabTransEstimator.Controls.Add(this.groupBox4);
            this.tabTransEstimator.Controls.Add(this.groupBox3);
            this.tabTransEstimator.Controls.Add(this.groupBox2);
            this.tabTransEstimator.Location = new System.Drawing.Point(4, 22);
            this.tabTransEstimator.Name = "tabTransEstimator";
            this.tabTransEstimator.Padding = new System.Windows.Forms.Padding(3);
            this.tabTransEstimator.Size = new System.Drawing.Size(624, 550);
            this.tabTransEstimator.TabIndex = 1;
            this.tabTransEstimator.Text = "Transformation estimation settings";
            this.tabTransEstimator.UseVisualStyleBackColor = true;
            // 
            // groupBox4
            // 
            this.groupBox4.Controls.Add(this.label4);
            this.groupBox4.Controls.Add(this.textBoxCorrQuality);
            this.groupBox4.Location = new System.Drawing.Point(8, 356);
            this.groupBox4.Name = "groupBox4";
            this.groupBox4.Size = new System.Drawing.Size(608, 72);
            this.groupBox4.TabIndex = 21;
            this.groupBox4.TabStop = false;
            this.groupBox4.Text = "Correspondences";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(16, 32);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(131, 13);
            this.label4.TabIndex = 13;
            this.label4.Text = "&Correspondence condition";
            // 
            // textBoxCorrQuality
            // 
            this.textBoxCorrQuality.Location = new System.Drawing.Point(176, 29);
            this.textBoxCorrQuality.Name = "textBoxCorrQuality";
            this.textBoxCorrQuality.Size = new System.Drawing.Size(80, 20);
            this.textBoxCorrQuality.TabIndex = 12;
            this.textBoxCorrQuality.Text = "0.75";
            //this.textBoxCorrQuality.TextChanged += new System.EventHandler(this.textBoxCorrQuality_TextChanged);
            // 
            // groupBox3
            // 
            this.groupBox3.Controls.Add(this.labelBlobScales);
            this.groupBox3.Controls.Add(this.textBoxBlobScales);
            this.groupBox3.Controls.Add(this.label6);
            this.groupBox3.Controls.Add(this.comboBoxPatchSize);
            this.groupBox3.Controls.Add(this.label5);
            this.groupBox3.Controls.Add(this.comboBoxDescriptor);
            this.groupBox3.Controls.Add(this.label3);
            this.groupBox3.Controls.Add(this.textBoxNumCorners);
            this.groupBox3.Controls.Add(this.label2);
            this.groupBox3.Controls.Add(this.comboBoxFeatureDetector);
            this.groupBox3.Location = new System.Drawing.Point(8, 150);
            this.groupBox3.Name = "groupBox3";
            this.groupBox3.Size = new System.Drawing.Size(608, 200);
            this.groupBox3.TabIndex = 20;
            this.groupBox3.TabStop = false;
            this.groupBox3.Text = "Features and Descriptors";
            // 
            // labelBlobScales
            // 
            this.labelBlobScales.AutoSize = true;
            this.labelBlobScales.Location = new System.Drawing.Point(16, 96);
            this.labelBlobScales.Name = "labelBlobScales";
            this.labelBlobScales.Size = new System.Drawing.Size(144, 13);
            this.labelBlobScales.TabIndex = 21;
            this.labelBlobScales.Text = "Num of &blob detection scales";
            // 
            // textBoxBlobScales
            // 
            this.textBoxBlobScales.Location = new System.Drawing.Point(176, 96);
            this.textBoxBlobScales.Name = "textBoxBlobScales";
            this.textBoxBlobScales.Size = new System.Drawing.Size(104, 20);
            this.textBoxBlobScales.TabIndex = 13;
            this.textBoxBlobScales.Text = "2";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(16, 160);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(56, 13);
            this.label6.TabIndex = 19;
            this.label6.Text = "&Patch size";
            // 
            // comboBoxPatchSize
            // 
            this.comboBoxPatchSize.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxPatchSize.FormattingEnabled = true;
            this.comboBoxPatchSize.Items.AddRange(new object[] {
            "9x9",
            "11x11",
            "13x13",
            "15x15"});
            this.comboBoxPatchSize.Location = new System.Drawing.Point(176, 160);
            this.comboBoxPatchSize.Name = "comboBoxPatchSize";
            this.comboBoxPatchSize.Size = new System.Drawing.Size(104, 21);
            this.comboBoxPatchSize.TabIndex = 18;
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(16, 128);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(107, 13);
            this.label5.TabIndex = 17;
            this.label5.Text = "Patch &descriptor type";
            // 
            // comboBoxDescriptor
            // 
            this.comboBoxDescriptor.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxDescriptor.FormattingEnabled = true;
            this.comboBoxDescriptor.Items.AddRange(new object[] {
            "Patches",
            "Oriented patches",
            "Normalised patches",
            "Oriented normalised patches"});
            this.comboBoxDescriptor.Location = new System.Drawing.Point(176, 128);
            this.comboBoxDescriptor.Name = "comboBoxDescriptor";
            this.comboBoxDescriptor.Size = new System.Drawing.Size(288, 21);
            this.comboBoxDescriptor.TabIndex = 16;
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(16, 64);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(129, 13);
            this.label3.TabIndex = 13;
            this.label3.Text = "&Target number of features";
            // 
            // textBoxNumCorners
            // 
            this.textBoxNumCorners.Location = new System.Drawing.Point(176, 64);
            this.textBoxNumCorners.Name = "textBoxNumCorners";
            this.textBoxNumCorners.Size = new System.Drawing.Size(104, 20);
            this.textBoxNumCorners.TabIndex = 12;
            this.textBoxNumCorners.Text = "150";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(16, 32);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(117, 13);
            this.label2.TabIndex = 11;
            this.label2.Text = "&Salient feature detector";
            // 
            // comboBoxFeatureDetector
            // 
            this.comboBoxFeatureDetector.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxFeatureDetector.FormattingEnabled = true;
            this.comboBoxFeatureDetector.Items.AddRange(new object[] {
            "Shi-Tomasi corners",
            "Shi-Tomasi corners with subpixel refinement",
            "Difference-of-Gaussian Blobs"});
            this.comboBoxFeatureDetector.Location = new System.Drawing.Point(176, 32);
            this.comboBoxFeatureDetector.Name = "comboBoxFeatureDetector";
            this.comboBoxFeatureDetector.Size = new System.Drawing.Size(288, 21);
            this.comboBoxFeatureDetector.TabIndex = 10;
            this.comboBoxFeatureDetector.SelectedIndexChanged += new System.EventHandler(this.comboBoxFeatureDetector_SelectedIndexChanged);
            // 
            // groupBox2
            // 
            this.groupBox2.Controls.Add(this.labelLMiters);
            this.groupBox2.Controls.Add(this.textBoxLMIterations);
            this.groupBox2.Controls.Add(this.checkBoxLM);
            this.groupBox2.Controls.Add(this.comboBoxTransType);
            this.groupBox2.Controls.Add(this.label1);
            this.groupBox2.Location = new System.Drawing.Point(8, 16);
            this.groupBox2.Name = "groupBox2";
            this.groupBox2.Size = new System.Drawing.Size(608, 128);
            this.groupBox2.TabIndex = 19;
            this.groupBox2.TabStop = false;
            this.groupBox2.Text = "Transformation";
            // 
            // labelLMiters
            // 
            this.labelLMiters.AutoSize = true;
            this.labelLMiters.Location = new System.Drawing.Point(16, 96);
            this.labelLMiters.Name = "labelLMiters";
            this.labelLMiters.Size = new System.Drawing.Size(122, 13);
            this.labelLMiters.TabIndex = 21;
            this.labelLMiters.Text = "Max &number of iterations";
            // 
            // textBoxLMIterations
            // 
            this.textBoxLMIterations.Location = new System.Drawing.Point(176, 96);
            this.textBoxLMIterations.Name = "textBoxLMIterations";
            this.textBoxLMIterations.Size = new System.Drawing.Size(80, 20);
            this.textBoxLMIterations.TabIndex = 20;
            this.textBoxLMIterations.Text = "8";
            // 
            // checkBoxLM
            // 
            this.checkBoxLM.AutoSize = true;
            this.checkBoxLM.Checked = true;
            this.checkBoxLM.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBoxLM.Location = new System.Drawing.Point(16, 64);
            this.checkBoxLM.Name = "checkBoxLM";
            this.checkBoxLM.Size = new System.Drawing.Size(309, 17);
            this.checkBoxLM.TabIndex = 19;
            this.checkBoxLM.Text = "Perform transformation &refinement (partial bundle adjustment)";
            this.checkBoxLM.UseVisualStyleBackColor = true;
            this.checkBoxLM.CheckedChanged += new System.EventHandler(this.checkBoxLM_CheckedChanged);
            // 
            // comboBoxTransType
            // 
            this.comboBoxTransType.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxTransType.FormattingEnabled = true;
            this.comboBoxTransType.Items.AddRange(new object[] {
            "Perspective transform",
            "Affine transform",
            "Similarity transform"});
            this.comboBoxTransType.Location = new System.Drawing.Point(176, 32);
            this.comboBoxTransType.Name = "comboBoxTransType";
            this.comboBoxTransType.Size = new System.Drawing.Size(288, 21);
            this.comboBoxTransType.TabIndex = 4;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(16, 32);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(100, 13);
            this.label1.TabIndex = 5;
            this.label1.Text = "&Transformation type";
            // 
            // tabRenderer
            // 
            this.tabRenderer.Controls.Add(this.groupBox5);
            this.tabRenderer.Location = new System.Drawing.Point(4, 22);
            this.tabRenderer.Name = "tabRenderer";
            this.tabRenderer.Size = new System.Drawing.Size(624, 550);
            this.tabRenderer.TabIndex = 2;
            this.tabRenderer.Text = "Renderer settings";
            this.tabRenderer.UseVisualStyleBackColor = true;
            // 
            // groupBox5
            // 
            this.groupBox5.Controls.Add(this.label13);
            this.groupBox5.Controls.Add(this.label12);
            this.groupBox5.Controls.Add(this.textBoxFullFrameUpdateFreq);
            this.groupBox5.Controls.Add(this.checkBoxSkipFrames);
            this.groupBox5.Controls.Add(this.label7);
            this.groupBox5.Controls.Add(this.comboBoxWarpMethod);
            this.groupBox5.Controls.Add(this.labelFeatherRad);
            this.groupBox5.Controls.Add(this.textBoxFeatherRad);
            this.groupBox5.Controls.Add(this.labelDijkstraScale);
            this.groupBox5.Controls.Add(this.textBoxDijkstraScale);
            this.groupBox5.Controls.Add(this.label11);
            this.groupBox5.Controls.Add(this.textBoxMosaicY);
            this.groupBox5.Controls.Add(this.label10);
            this.groupBox5.Controls.Add(this.textBoxMosaicX);
            this.groupBox5.Controls.Add(this.label9);
            this.groupBox5.Controls.Add(this.textBoxMaxFrames);
            this.groupBox5.Controls.Add(this.checkBoxIncrementalRenderer);
            this.groupBox5.Controls.Add(this.label8);
            this.groupBox5.Controls.Add(this.comboBoxRenderer);
            this.groupBox5.Location = new System.Drawing.Point(8, 16);
            this.groupBox5.Name = "groupBox5";
            this.groupBox5.Size = new System.Drawing.Size(608, 436);
            this.groupBox5.TabIndex = 22;
            this.groupBox5.TabStop = false;
            this.groupBox5.Text = "Renderer";
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(245, 164);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(38, 13);
            this.label13.TabIndex = 38;
            this.label13.Text = "frames";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(32, 162);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(138, 13);
            this.label12.TabIndex = 37;
            this.label12.Text = "Perform full rendering every ";
            // 
            // textBoxFullFrameUpdateFreq
            // 
            this.textBoxFullFrameUpdateFreq.Location = new System.Drawing.Point(176, 159);
            this.textBoxFullFrameUpdateFreq.Name = "textBoxFullFrameUpdateFreq";
            this.textBoxFullFrameUpdateFreq.Size = new System.Drawing.Size(64, 20);
            this.textBoxFullFrameUpdateFreq.TabIndex = 36;
            this.textBoxFullFrameUpdateFreq.Text = "7";
            // 
            // checkBoxSkipFrames
            // 
            this.checkBoxSkipFrames.AutoSize = true;
            this.checkBoxSkipFrames.Location = new System.Drawing.Point(16, 202);
            this.checkBoxSkipFrames.Name = "checkBoxSkipFrames";
            this.checkBoxSkipFrames.Size = new System.Drawing.Size(208, 17);
            this.checkBoxSkipFrames.TabIndex = 35;
            this.checkBoxSkipFrames.Text = "Use &selective rendering (larger mosaic)";
            this.checkBoxSkipFrames.UseVisualStyleBackColor = true;
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(16, 242);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(131, 13);
            this.label7.TabIndex = 34;
            this.label7.Text = "&Warp interpolation method";
            // 
            // comboBoxWarpMethod
            // 
            this.comboBoxWarpMethod.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxWarpMethod.FormattingEnabled = true;
            this.comboBoxWarpMethod.Items.AddRange(new object[] {
            "Nearest neighbour (fastest)",
            "Bilinear (best)"});
            this.comboBoxWarpMethod.Location = new System.Drawing.Point(176, 242);
            this.comboBoxWarpMethod.Name = "comboBoxWarpMethod";
            this.comboBoxWarpMethod.Size = new System.Drawing.Size(288, 21);
            this.comboBoxWarpMethod.TabIndex = 33;
            // 
            // labelFeatherRad
            // 
            this.labelFeatherRad.AutoSize = true;
            this.labelFeatherRad.Location = new System.Drawing.Point(16, 96);
            this.labelFeatherRad.Name = "labelFeatherRad";
            this.labelFeatherRad.Size = new System.Drawing.Size(74, 13);
            this.labelFeatherRad.TabIndex = 25;
            this.labelFeatherRad.Text = "&Feather radius";
            // 
            // textBoxFeatherRad
            // 
            this.textBoxFeatherRad.Location = new System.Drawing.Point(176, 96);
            this.textBoxFeatherRad.Name = "textBoxFeatherRad";
            this.textBoxFeatherRad.Size = new System.Drawing.Size(64, 20);
            this.textBoxFeatherRad.TabIndex = 26;
            this.textBoxFeatherRad.Text = "16";
            // 
            // labelDijkstraScale
            // 
            this.labelDijkstraScale.AutoSize = true;
            this.labelDijkstraScale.Location = new System.Drawing.Point(16, 64);
            this.labelDijkstraScale.Name = "labelDijkstraScale";
            this.labelDijkstraScale.Size = new System.Drawing.Size(120, 13);
            this.labelDijkstraScale.TabIndex = 23;
            this.labelDijkstraScale.Text = "&Dijkstra scale parameter";
            // 
            // textBoxDijkstraScale
            // 
            this.textBoxDijkstraScale.Location = new System.Drawing.Point(176, 64);
            this.textBoxDijkstraScale.Name = "textBoxDijkstraScale";
            this.textBoxDijkstraScale.Size = new System.Drawing.Size(64, 20);
            this.textBoxDijkstraScale.TabIndex = 24;
            this.textBoxDijkstraScale.Text = "4";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Font = new System.Drawing.Font("Microsoft Sans Serif", 11.25F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.label11.Location = new System.Drawing.Point(248, 290);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(16, 18);
            this.label11.TabIndex = 30;
            this.label11.Text = "x";
            // 
            // textBoxMosaicY
            // 
            this.textBoxMosaicY.Location = new System.Drawing.Point(272, 290);
            this.textBoxMosaicY.Name = "textBoxMosaicY";
            this.textBoxMosaicY.Size = new System.Drawing.Size(64, 20);
            this.textBoxMosaicY.TabIndex = 30;
            this.textBoxMosaicY.Text = "1000";
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(16, 290);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(96, 13);
            this.label10.TabIndex = 28;
            this.label10.Text = "&Mosaic dimensions";
            // 
            // textBoxMosaicX
            // 
            this.textBoxMosaicX.Location = new System.Drawing.Point(176, 290);
            this.textBoxMosaicX.Name = "textBoxMosaicX";
            this.textBoxMosaicX.Size = new System.Drawing.Size(64, 20);
            this.textBoxMosaicX.TabIndex = 29;
            this.textBoxMosaicX.Text = "1000";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(16, 320);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(115, 13);
            this.label9.TabIndex = 31;
            this.label9.Text = "Ma&x frames per mosaic";
            // 
            // textBoxMaxFrames
            // 
            this.textBoxMaxFrames.Location = new System.Drawing.Point(176, 320);
            this.textBoxMaxFrames.Name = "textBoxMaxFrames";
            this.textBoxMaxFrames.Size = new System.Drawing.Size(64, 20);
            this.textBoxMaxFrames.TabIndex = 32;
            this.textBoxMaxFrames.Text = "30";
            // 
            // checkBoxIncrementalRenderer
            // 
            this.checkBoxIncrementalRenderer.AutoSize = true;
            this.checkBoxIncrementalRenderer.Checked = true;
            this.checkBoxIncrementalRenderer.CheckState = System.Windows.Forms.CheckState.Checked;
            this.checkBoxIncrementalRenderer.Location = new System.Drawing.Point(16, 136);
            this.checkBoxIncrementalRenderer.Name = "checkBoxIncrementalRenderer";
            this.checkBoxIncrementalRenderer.Size = new System.Drawing.Size(184, 17);
            this.checkBoxIncrementalRenderer.TabIndex = 27;
            this.checkBoxIncrementalRenderer.Text = "Use &incremental rendering (faster)";
            this.checkBoxIncrementalRenderer.UseVisualStyleBackColor = true;
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(16, 32);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(51, 13);
            this.label8.TabIndex = 23;
            this.label8.Text = "&Renderer";
            // 
            // comboBoxRenderer
            // 
            this.comboBoxRenderer.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxRenderer.FormattingEnabled = true;
            this.comboBoxRenderer.Items.AddRange(new object[] {
            "Basic Renderer",
            "Dijkstra Shortest Path Cut Renderer",
            "Feathered Renderer",
            "Multi-Scale Feathered Renderer"});
            this.comboBoxRenderer.Location = new System.Drawing.Point(176, 32);
            this.comboBoxRenderer.Name = "comboBoxRenderer";
            this.comboBoxRenderer.Size = new System.Drawing.Size(288, 21);
            this.comboBoxRenderer.TabIndex = 22;
            this.comboBoxRenderer.SelectedIndexChanged += new System.EventHandler(this.comboBoxRenderer_SelectedIndexChanged);
            // 
            // tabPageEval
            // 
            this.tabPageEval.Controls.Add(this.groupBox6);
            this.tabPageEval.Location = new System.Drawing.Point(4, 22);
            this.tabPageEval.Name = "tabPageEval";
            this.tabPageEval.Size = new System.Drawing.Size(624, 550);
            this.tabPageEval.TabIndex = 3;
            this.tabPageEval.Text = "Evaluation method";
            this.tabPageEval.UseVisualStyleBackColor = true;
            // 
            // groupBox6
            // 
            this.groupBox6.Controls.Add(this.label18);
            this.groupBox6.Controls.Add(this.comboBoxEvalFunction);
            this.groupBox6.Location = new System.Drawing.Point(8, 16);
            this.groupBox6.Name = "groupBox6";
            this.groupBox6.Size = new System.Drawing.Size(608, 72);
            this.groupBox6.TabIndex = 23;
            this.groupBox6.TabStop = false;
            this.groupBox6.Text = "Mosaic quality evaluation function";
            // 
            // label18
            // 
            this.label18.AutoSize = true;
            this.label18.Location = new System.Drawing.Point(16, 32);
            this.label18.Name = "label18";
            this.label18.Size = new System.Drawing.Size(98, 13);
            this.label18.TabIndex = 23;
            this.label18.Text = "&Evaluation function";
            // 
            // comboBoxEvalFunction
            // 
            this.comboBoxEvalFunction.DropDownStyle = System.Windows.Forms.ComboBoxStyle.DropDownList;
            this.comboBoxEvalFunction.FormattingEnabled = true;
            this.comboBoxEvalFunction.Items.AddRange(new object[] {
            "None",
            "Sum-of-Squared-Differences evaluation function"});
            this.comboBoxEvalFunction.Location = new System.Drawing.Point(176, 32);
            this.comboBoxEvalFunction.Name = "comboBoxEvalFunction";
            this.comboBoxEvalFunction.Size = new System.Drawing.Size(288, 21);
            this.comboBoxEvalFunction.TabIndex = 22;
            // 
            // groupBox7
            // 
            this.groupBox7.Controls.Add(this.label15);
            this.groupBox7.Controls.Add(this.label14);
            this.groupBox7.Controls.Add(this.textBoxMaxSearchForTransform);
            this.groupBox7.Controls.Add(this.textBoxMaxRansacIters);
            this.groupBox7.Location = new System.Drawing.Point(9, 441);
            this.groupBox7.Name = "groupBox7";
            this.groupBox7.Size = new System.Drawing.Size(606, 99);
            this.groupBox7.TabIndex = 22;
            this.groupBox7.TabStop = false;
            this.groupBox7.Text = "Performance";
            // 
            // textBoxMaxRansacIters
            // 
            this.textBoxMaxRansacIters.Location = new System.Drawing.Point(175, 29);
            this.textBoxMaxRansacIters.Name = "textBoxMaxRansacIters";
            this.textBoxMaxRansacIters.Size = new System.Drawing.Size(79, 20);
            this.textBoxMaxRansacIters.TabIndex = 0;
            this.textBoxMaxRansacIters.Text = "1000";
            // 
            // textBoxMaxSearchForTransform
            // 
            this.textBoxMaxSearchForTransform.Location = new System.Drawing.Point(176, 55);
            this.textBoxMaxSearchForTransform.Name = "textBoxMaxSearchForTransform";
            this.textBoxMaxSearchForTransform.Size = new System.Drawing.Size(79, 20);
            this.textBoxMaxSearchForTransform.TabIndex = 1;
            this.textBoxMaxSearchForTransform.Text = "10";
            //this.textBoxMaxSearchForTransform.TextChanged += new System.EventHandler(this.textBoxMaxSearchForTransform_TextChanged);
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.Location = new System.Drawing.Point(15, 32);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(119, 13);
            this.label14.TabIndex = 2;
            this.label14.Text = "Max RANSAC iterations";
            // 
            // label15
            // 
            this.label15.AutoSize = true;
            this.label15.Location = new System.Drawing.Point(15, 58);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(108, 13);
            this.label15.TabIndex = 3;
            this.label15.Text = "Max frames to search";
            // 
            // SettingsGUI
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(647, 630);
            this.Controls.Add(this.tabControlMain);
            this.Controls.Add(this.buttonExit);
            this.Controls.Add(this.buttonGo);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedDialog;
            this.MaximizeBox = false;
            this.Name = "SettingsGUI";
            this.Text = "DTA Video Mosaic";
            this.tabControlMain.ResumeLayout(false);
            this.tabVidSource.ResumeLayout(false);
            this.groupBoxOutput.ResumeLayout(false);
            this.groupBoxOutput.PerformLayout();
            this.groupBox1.ResumeLayout(false);
            this.groupBox1.PerformLayout();
            this.tabTransEstimator.ResumeLayout(false);
            this.groupBox4.ResumeLayout(false);
            this.groupBox4.PerformLayout();
            this.groupBox3.ResumeLayout(false);
            this.groupBox3.PerformLayout();
            this.groupBox2.ResumeLayout(false);
            this.groupBox2.PerformLayout();
            this.tabRenderer.ResumeLayout(false);
            this.groupBox5.ResumeLayout(false);
            this.groupBox5.PerformLayout();
            this.tabPageEval.ResumeLayout(false);
            this.groupBox6.ResumeLayout(false);
            this.groupBox6.PerformLayout();
            this.groupBox7.ResumeLayout(false);
            this.groupBox7.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button buttonGo;
        private System.Windows.Forms.Button buttonExit;
        private System.Windows.Forms.TabControl tabControlMain;
        private System.Windows.Forms.TabPage tabVidSource;
        private System.Windows.Forms.TabPage tabTransEstimator;
        private System.Windows.Forms.TabPage tabRenderer;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.RadioButton radioButtonDir;
        private System.Windows.Forms.RadioButton radioButtonVideo;
        private System.Windows.Forms.RadioButton radioButtonCam;
        private System.Windows.Forms.TextBox textBoxDir;
        private System.Windows.Forms.TextBox textBoxVidFile;
        private System.Windows.Forms.Button buttonBrowseVidFile;
        private System.Windows.Forms.Button buttonBrowseDir;
        private System.Windows.Forms.GroupBox groupBoxOutput;
        private System.Windows.Forms.Button buttonBrowseSaveIm;
        private System.Windows.Forms.TextBox textBoxSaveMosaic;
        private System.Windows.Forms.RadioButton radioButtonSaveMosaic;
        private System.Windows.Forms.RadioButton radioButtonDisplayOnly;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.ComboBox comboBoxTransType;
        private System.Windows.Forms.GroupBox groupBox2;
        private System.Windows.Forms.Label labelLMiters;
        private System.Windows.Forms.TextBox textBoxLMIterations;
        private System.Windows.Forms.CheckBox checkBoxLM;
        private System.Windows.Forms.GroupBox groupBox3;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.ComboBox comboBoxPatchSize;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.ComboBox comboBoxDescriptor;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox textBoxNumCorners;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.ComboBox comboBoxFeatureDetector;
        private System.Windows.Forms.GroupBox groupBox4;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.TextBox textBoxCorrQuality;
        private System.Windows.Forms.GroupBox groupBox5;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.TextBox textBoxMaxFrames;
        private System.Windows.Forms.CheckBox checkBoxIncrementalRenderer;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.ComboBox comboBoxRenderer;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.TextBox textBoxMosaicX;
        private System.Windows.Forms.TextBox textBoxMosaicY;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label labelBlobScales;
        private System.Windows.Forms.TextBox textBoxBlobScales;
        private System.Windows.Forms.TabPage tabPageEval;
        private System.Windows.Forms.Label labelFeatherRad;
        private System.Windows.Forms.TextBox textBoxFeatherRad;
        private System.Windows.Forms.Label labelDijkstraScale;
        private System.Windows.Forms.TextBox textBoxDijkstraScale;
        private System.Windows.Forms.GroupBox groupBox6;
        private System.Windows.Forms.Label label18;
        private System.Windows.Forms.ComboBox comboBoxEvalFunction;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.ComboBox comboBoxWarpMethod;
        private System.Windows.Forms.CheckBox checkBoxSkipFrames;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.TextBox textBoxFullFrameUpdateFreq;
        private System.Windows.Forms.GroupBox groupBox7;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.TextBox textBoxMaxSearchForTransform;
        private System.Windows.Forms.TextBox textBoxMaxRansacIters;
    }
}

