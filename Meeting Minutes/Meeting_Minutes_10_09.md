# Team 9 Meeting Minutes 
# Meeting Details:
## Date: 
October 9th, 2024

## Location:
Zoom: https://ubc.zoom.us/j/67359231965?pwd=1QOWWuFTTjXQcla2T32GQ5MoyAUHhG.1

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- Research question: Is the microbiome a mediator of the increased prevalence of asthma in patients with MS compared to the general public?
    - More specific question: Determine whether the microbiota differs in RRMS and PMS patients with and without asthma.
- Metadata discussion
    - Variables retained: sample-id, site_x, sex, age, year_of_onset, disease, disease_course, disease_duration, diet_no_special_needs,         diet_special_needs, birth_method, breastfeeding, allergies, allergy_specific, asthma, eczema, ms_family, nsaids, nsaids_specifics,        rxmeds, rxmeds_number, otc_meds, otc_number, smoke
    - Are there any variables missing that would be valuable for our analysis?
          - We haven't included MS treatment variables
    - Are there any variables retained that are not important for the analysis?
          - We wanted to keep variables that would have potential impacts on the microbiome.
    - Concerns: how do we control for confounding variables?
    - There aren't duplicates for every sample - should we filter out the duplicates?
 - MS + asthma prevalence
    - 76% RRMS (437/576)
    - PMS: 12% SPMS (68/576) + 12% PPMS (71/576)
    - 2 patients with PPMS have asthma
    - 38 patients with RRMS have asthma
    - Should we filter out patients with PMS?
- Sample sizes for each group:
    - RRMS + Asthma: 38
    - RRMS + no asthma: 38
    - no MS + no asthma: 38
    - We haven't randomly sampled yet
- Make the manifest file - TA assistance? 
- Filter the metadata in R once variables are confirmed


