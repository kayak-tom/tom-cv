using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
//using grc;

namespace DTA_Mosaicing_GUI
{
    public partial class SettingsGUI : Form
    {
        public SettingsGUI()
        {
            InitializeComponent();

            comboBoxDescriptor.SelectedIndex = ((int)grc.Enums.eDescriptorInvarianceType.eOrientedNormalised);
            comboBoxFeatureDetector.SelectedIndex = ((int)grc.Enums.eFeatureDetector.eShiTomasiSubpix);
            comboBoxPatchSize.SelectedIndex = ((int)grc.Enums.eDescriptorSize.e13x13);
            comboBoxRenderer.SelectedIndex = ((int)grc.Enums.eRenderer.eDijkstraRenderer);
            comboBoxTransType.SelectedIndex = ((int)grc.Enums.eTransformType.ePerspective);
            comboBoxEvalFunction.SelectedIndex = ((int)grc.Enums.eRendererEvalFunction.eNoEvaluation);
            comboBoxWarpMethod.SelectedIndex = ((int)grc.Enums.eWarpMethod.eWarpNN);

            checkBoxLM_CheckedChanged(null, null);
            doFileboxEnable();
        }

        private bool validateIOSettings()
        {
            tabControlMain.SelectedIndex = 0;

            if (radioButtonVideo.Checked)
                if (!System.IO.File.Exists(textBoxVidFile.Text))
                {
                    textBoxVidFile.SelectAll();
                    textBoxVidFile.Focus();
                    MessageBox.Show("Please select a video file", "File not found");
                    return false;
                }
            if (radioButtonDir.Checked)
                if (!System.IO.Directory.Exists(textBoxDir.Text))
                {
                    textBoxDir.SelectAll();
                    textBoxDir.Focus();
                    MessageBox.Show("Please select a folder containing image files", "Folder not found");
                    return false;
                }
            if (radioButtonSaveMosaic.Checked)
                if (!System.IO.Directory.Exists(textBoxSaveMosaic.Text))
                {
                    textBoxSaveMosaic.SelectAll();
                    textBoxSaveMosaic.Focus();
                    MessageBox.Show("Please select a folder to save mosaics to", "Folder not found");
                    return false;
                }
            return true;
        }
        private bool validateTransformSettings()
        {
            tabControlMain.SelectedIndex = 1;

            int numCorners;
            if (!int.TryParse(textBoxNumCorners.Text, out numCorners)
                || numCorners < 5
                || numCorners > 10000 )
            {
                textBoxNumCorners.SelectAll();
                textBoxNumCorners.Focus();
                MessageBox.Show("Target number of corners must be an integer in the range 5-10000");
                return false;
            }

            double corrCondition;
            if (!double.TryParse(textBoxCorrQuality.Text, out corrCondition)
                || corrCondition < 0.1
                || corrCondition > 1.0 )
            {
                textBoxCorrQuality.SelectAll();
                textBoxCorrQuality.Focus();
                MessageBox.Show("Correspondence condition must be between 0.1 (only very strong matches accepted) and 1.0 (all correspondences accepted)");
                return false;
            }

            int maxRansacIterations;
            if (!int.TryParse(textBoxMaxRansacIters.Text, out maxRansacIterations)
                || maxRansacIterations < 100
                || maxRansacIterations > 10000)
            {
                textBoxMaxRansacIters.SelectAll();
                textBoxMaxRansacIters.Focus();
                MessageBox.Show("Maximum ransac iterations must be an integer in the range 100-10000");
                return false;
            }

            int maxSearchForTransform;
            if (!int.TryParse(textBoxMaxSearchForTransform.Text, out maxSearchForTransform)
                || maxSearchForTransform < 1
                || maxSearchForTransform > 50)
            {
                textBoxMaxSearchForTransform.SelectAll();
                textBoxMaxSearchForTransform.Focus();
                MessageBox.Show("Maximum number of frames to search for transform must be an integer in the range 1-50");
                return false;
            }

            if (checkBoxLM.Checked)
            {
                int numLMIters;
                if (!int.TryParse(textBoxLMIterations.Text, out numLMIters)
                    || numLMIters < 1
                    || numLMIters > 1000)
                {
                    textBoxLMIterations.SelectAll();
                    textBoxLMIterations.Focus();
                    MessageBox.Show("Number of Levenberg-Marquardt iterations must be an integer in the range 1-1000");
                    return false;
                }
            }

            if (comboBoxFeatureDetector.SelectedIndex == (int)grc.Enums.eFeatureDetector.eDoGBlob)
            {
                int numBlobScales;
                if (!int.TryParse(textBoxBlobScales.Text, out numBlobScales)
                    || numBlobScales < 1
                    || numBlobScales > 4)
                {
                    textBoxBlobScales.SelectAll();
                    textBoxBlobScales.Focus();
                    MessageBox.Show("Number of blob-detection scales must be an integer in the range 1-4");
                    return false;
                }
            }
            
            return true;
        }
        private bool validateRendererSettings()
        {
            tabControlMain.SelectedIndex = 2;
            if (comboBoxFeatureDetector.SelectedIndex == (int)grc.Enums.eRenderer.eFeatheredRenderer || comboBoxFeatureDetector.SelectedIndex == (int)grc.Enums.eRenderer.eMultiScaleRenderer)
            {
                int featherRad;
                if (!int.TryParse(textBoxFeatherRad.Text, out featherRad)
                    || featherRad < 1
                    || featherRad > 64)
                {
                    textBoxFeatherRad.SelectAll();
                    textBoxFeatherRad.Focus();
                    MessageBox.Show("Feathering radius must be an integer in the range 1-64");
                    return false;
                }
            }
            
            if (comboBoxFeatureDetector.SelectedIndex == (int)grc.Enums.eRenderer.eDijkstraRenderer)
            {
                int dijkstraScale;
                if (!int.TryParse(textBoxDijkstraScale.Text, out dijkstraScale)
                    || dijkstraScale < 1
                    || dijkstraScale > 8)
                {
                    textBoxDijkstraScale.SelectAll();
                    textBoxDijkstraScale.Focus();
                    MessageBox.Show("Dijkstra-cut scale must be an integer in the range 1-8");
                    return false;
                }
            }
            int maxFrames;
            if (!int.TryParse(textBoxMaxFrames.Text, out maxFrames)
                || maxFrames < 1
                || maxFrames > 1000)
            {
                textBoxMaxFrames.SelectAll();
                textBoxMaxFrames.Focus();
                MessageBox.Show("The maximum number of frames to render must be an integer in the range 1-1000");
                return false;
            }
            int mosaicX, mosaicY;
            if (!int.TryParse(textBoxMosaicX.Text, out mosaicX)
                || mosaicX < 100
                || mosaicX > 10000
                || !int.TryParse(textBoxMosaicY.Text, out mosaicY)
                || mosaicY < 100
                || mosaicY > 10000)
            {
                textBoxMosaicX.SelectAll();
                textBoxMosaicX.Focus();
                MessageBox.Show("Each mosaic dimention must be an integer in the range 100-10000");
                return false;
            }

            int fullRenderFreq;
            
            if (!int.TryParse(textBoxFullFrameUpdateFreq.Text, out fullRenderFreq)
                || fullRenderFreq < 3
                || fullRenderFreq > 30)
            {
                textBoxFullFrameUpdateFreq.SelectAll();
                textBoxFullFrameUpdateFreq.Focus();
                MessageBox.Show("Frame update frequency must be an integer in the range 3-30");
                return false;
            }
            
            

            return true;
        }
        private bool validateForm()
        {
            int tabSel = tabControlMain.SelectedIndex;
            bool bSuccess = validateIOSettings() && validateTransformSettings() && validateRendererSettings();
            if (bSuccess)
                tabControlMain.SelectedIndex = tabSel;
            return bSuccess;
        }

        class configFile : System.IO.StreamWriter
        {
            public configFile(string cfgPath) : base(cfgPath) {}
            
            public void addParam(string name, object val)
            {
                string value = val.ToString();
                base.WriteLine(name + '=' + value);
            }
        }

        int getVidSource()
        {
            if (radioButtonCam.Checked)
                return (int)grc.Enums.eVideoSource.eAttachedCam;
            if (radioButtonDir.Checked)
                return (int)grc.Enums.eVideoSource.eImageDirectory;
            if (radioButtonVideo.Checked)
                return (int)grc.Enums.eVideoSource.eVideoFile;

            throw new ApplicationException();
        }
        int getMosaicDestination()
        {
            if (radioButtonSaveMosaic.Checked)
                return (int)grc.Enums.eDisplayMode.eSaveImages;
            if (radioButtonDisplayOnly.Checked)
                return (int)grc.Enums.eDisplayMode.eDisplayOnly;

            throw new ApplicationException();
        }

        /*private void initComboBoxes()
        {
            comboBoxTransType.Items
        }*/

        private void saveConfig(string cfgPath)
        {
//            System.IO.TextWriter config = new System.IO.StreamWriter(cfgPath);
            configFile config = new configFile(cfgPath);

            // I/O settings
            config.addParam("VideoSource", getVidSource());
            config.addParam("ImageDir", textBoxDir.Text );
            config.addParam("VideoFile", textBoxVidFile.Text);

            config.addParam("MosaicDestination", getMosaicDestination());
            config.addParam("MosaicSaveDir", textBoxSaveMosaic.Text);

            // Transform settings
            config.addParam("TransformType", comboBoxTransType.SelectedIndex);
            config.addParam("FeatureDetector", comboBoxFeatureDetector.SelectedIndex);
            config.addParam("NumFeatures", textBoxNumCorners.Text);
            config.addParam("NumBlobScales", textBoxBlobScales.Text);
            config.addParam("CorrespondenceQuality", textBoxCorrQuality.Text);
            config.addParam("LMIterations", textBoxLMIterations.Text);
            config.addParam("LM", checkBoxLM.Checked ? 1 : 0);
            config.addParam("PatchDescriptorType", comboBoxDescriptor.SelectedIndex);
            config.addParam("PatchDescriptorSize", comboBoxPatchSize.SelectedIndex);
            config.addParam("MaxRansacIters", textBoxMaxRansacIters.Text);
            config.addParam("MaxSearchForTransform", textBoxMaxSearchForTransform.Text);

            //Renderer settings
            config.addParam("MaxFrames", textBoxMaxFrames.Text);
            config.addParam("IncrementalRendering", checkBoxIncrementalRenderer.Checked ? 1 : 0);
            config.addParam("FullFrameUpdateFreq", textBoxFullFrameUpdateFreq.Text);
            config.addParam("SkipFrames", checkBoxSkipFrames.Checked ? 1 : 0);
            config.addParam("RendererType", comboBoxRenderer.SelectedIndex);
            config.addParam("FeatherRadius", textBoxFeatherRad.Text);
            config.addParam("DijkstraScale", textBoxDijkstraScale.Text);

            config.addParam("WarpMethod", comboBoxWarpMethod.SelectedIndex);

            config.addParam("MosaicX", textBoxMosaicX.Text);
            config.addParam("MosaicY", textBoxMosaicY.Text);
            
            //Evaluation settings
            config.addParam("EvaluationFunction", comboBoxEvalFunction.SelectedIndex);

            config.Close();
        }

        private void buttonExit_Click(object sender, EventArgs e)
        {
            //DTA_Mosaicing_GUI
            Application.Exit();
        }

        private void buttonGo_Click(object sender, EventArgs e)
        {
            try
            {

                if (!validateForm()) return;

                string path;
                path = System.IO.Path.GetDirectoryName(
                   System.Reflection.Assembly.GetExecutingAssembly().GetName().CodeBase);

                path = path.Replace('\\', '/');
                path = path.Replace("file:/", "");

                string cfgPath = path + "\\MosaicSettings.cfg";
                saveConfig(cfgPath);

                string mosaicerPath = path;
                if(path.ToLower().Contains("/release")) 
                    mosaicerPath = mosaicerPath + "/../../../Release" ;
                else if(path.ToLower().Contains("/debug")) 
                    mosaicerPath = mosaicerPath + "/../../../Debug" ;
                //Otherwise assume they're in the same folder
                mosaicerPath += "/DTA - Image Mosaicing 2.exe";

                if (!System.IO.File.Exists(mosaicerPath))
                {
                    MessageBox.Show("Can't find mosaicing application: Looking for \"" + mosaicerPath + "\"");
                    return;
                }

                System.Diagnostics.Process proc = new System.Diagnostics.Process();
                proc.EnableRaisingEvents = false;
                proc.StartInfo.FileName = mosaicerPath;// "C:/Documents and Settings/Tom Botterill/My Documents/GRC Work/dta-mosaicing/MosaicingSoftware/release/DTA - Image Mosaicing 2.exe"; // ;
                proc.StartInfo.Arguments = "\"" + cfgPath + "\"";
                proc.Start();
            }
            catch(ApplicationException apErr)
            {
                MessageBox.Show("Unknown application error: " + apErr.ToString());
            }
        }

        private void buttonBrowse_Click(object sender, EventArgs e)
        {
            OpenFileDialog fdlg = new OpenFileDialog();
            fdlg.InitialDirectory = textBoxVidFile.Text;
            fdlg.Filter = "All video files (*.avi;*.vfw;*.mpg)|*.avi;*.vfw;*.mpg|AVI Video (*.avi)|*.avi|VFW Video (*.vfw)|*.vfw|MPEG Video (*.mpg;*.mpeg)|*.mpg;*.mpeg|All files (*.*)|*.*";

            if (fdlg.ShowDialog() == DialogResult.OK)
            {
                textBoxVidFile.Text = fdlg.FileName;
            }
        }

        private void buttonBrowseDir_Click(object sender, EventArgs e)
        {
            FolderBrowserDialog fdlg = new FolderBrowserDialog();
            if (fdlg.ShowDialog() == DialogResult.OK)
            {
                textBoxDir.Text = fdlg.SelectedPath;
            }
        }

        private void buttonBrowseSaveIm_Click(object sender, EventArgs e)
        {
            FolderBrowserDialog fdlg = new FolderBrowserDialog();
            if (fdlg.ShowDialog() == DialogResult.OK)
            {
                textBoxSaveMosaic.Text = fdlg.SelectedPath;
            }
        }

        private void textBoxVidFile_TextChanged(object sender, EventArgs e)
        {
            radioButtonVideo.Checked = true;
        }

        private void textBoxDir_TextChanged(object sender, EventArgs e)
        {
            radioButtonDir.Checked = true;
        }

        private void textBoxSaveMosaic_TextChanged(object sender, EventArgs e)
        {
            radioButtonSaveMosaic.Checked = true;
        }

        private void checkBoxLM_CheckedChanged(object sender, EventArgs e)
        {
            textBoxLMIterations.Enabled = checkBoxLM.Checked;
            labelLMiters.Enabled = checkBoxLM.Checked;
        }

        private void comboBoxFeatureDetector_SelectedIndexChanged(object sender, EventArgs e)
        {
            bool enableBlob = (comboBoxFeatureDetector.SelectedIndex == (int)grc.Enums.eFeatureDetector.eDoGBlob);
            textBoxBlobScales.Enabled = enableBlob;
            labelBlobScales.Enabled = enableBlob;
        }

        private void comboBoxRenderer_SelectedIndexChanged(object sender, EventArgs e)
        {
            bool enableDijk = (comboBoxRenderer.SelectedIndex == (int)grc.Enums.eRenderer.eDijkstraRenderer);
            bool enableFeather = (comboBoxRenderer.SelectedIndex == (int)grc.Enums.eRenderer.eFeatheredRenderer)
                || (comboBoxRenderer.SelectedIndex == (int)grc.Enums.eRenderer.eMultiScaleRenderer);

            labelFeatherRad.Enabled = enableFeather;
            textBoxFeatherRad.Enabled = enableFeather;

            labelDijkstraScale.Enabled = enableDijk;
            textBoxDijkstraScale.Enabled = enableDijk;
        }
        private void doFileboxEnable()
        {
            bool enVidFile = radioButtonVideo.Checked;
            bool enImDir = radioButtonDir.Checked;
            bool enSaveDir = radioButtonSaveMosaic.Checked;

            textBoxSaveMosaic.Enabled = enSaveDir;
            buttonBrowseSaveIm.Enabled = enSaveDir;

            textBoxDir.Enabled = enImDir;
            buttonBrowseDir.Enabled = enImDir;

            textBoxVidFile.Enabled = enVidFile;
            buttonBrowseVidFile.Enabled = enVidFile;
        }
        private void radioButtonCam_CheckedChanged(object sender, EventArgs e)
        {
            doFileboxEnable();
        }

        private void radioButtonVideo_CheckedChanged(object sender, EventArgs e)
        {
            doFileboxEnable();
        }

        private void radioButtonDir_CheckedChanged(object sender, EventArgs e)
        {
            doFileboxEnable();
        }

        private void radioButtonDisplayOnly_CheckedChanged(object sender, EventArgs e)
        {
            doFileboxEnable();
        }

        private void radioButtonSaveMosaic_CheckedChanged(object sender, EventArgs e)
        {
            doFileboxEnable();
        }

    


    }
}