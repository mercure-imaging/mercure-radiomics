# **mercure-radiomics WORK IN PROGRESS**
<br>

Mercure module to run IBSI compliant radiomics feature extraction.

Mercure module to deploy the [LARS](https://github.com/haeggsti/lymphoma_classification_2023/) model for Lymphoma classification . This module runs as a docker container in mercure, it can be added to an existing mercure installation using docker tag : *mercureimaging/mercure-lars*. The module generates predictions of LARS-avg and LARS-max as described in the paper below.

*"Deep learning for 18F-fluorodeoxyglucose-PET-CT classification in patients with lymphoma: a dual-centre retrospective analysis", by I. Häggström et al., published in The Lancet Digital Health (2023).*

<br>

## **Installation**

### Add module to existing mercure installation
Follow instructions on [mercure website](https://mercure-imaging.org) on how to add a new module. Use the docker tag mercureimaging/mercure-lars.

<br>

### Build module for local testing, modification and development
1. Clone repo.
2. Add 'models' directory and add checkpoint files for top-10 trained models as found [here.](https://drive.google.com/drive/folders/1V-hhATi3zaqAiVyZ8_hgE3zhtSdt2HbV)
3. Add 'checkpoints' directory and add ResNet34 model file.
2. Build Docker container locally by running make (modify makefile with new docker tag as needed).
3. Test container :\
`docker run -it -v /input_data:/input -v /output_data:/output --env MERCURE_IN_DIR=/input  --env MERCURE_OUT_DIR=/output mercureimaging/mercure-lars`

<br>

## **Configuration**

The mercure-lars module requires little configuration in mercure. It runs on a single PET series and it is important that the processing rules are configured to filter out any surplus data received by mercure. More information on mercure rule configuration can be found [here.](https://mercure-imaging.org/docs/usage.html)

During preprocessing before inference, the module performs gaussian filtering by default ( set to 'true' ). This can be disabled in the settings via the Modules or Rules page as below: 
```
"gaussian_filter": "false" (default: true) 
```

<br>
<br>

## **Acknowledgments**

<br>
Models and inference code adapted from the lymphoma_classification_2023 repository: https://github.com/haeggsti/lymphoma_classification_2023
<br>



