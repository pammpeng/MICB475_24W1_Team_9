# Team 9 Meeting Minutes 
# Meeting Details:
## Date: 
October 9th, 2024

## Location:
Zoom: https://ubc.zoom.us/j/67359231965?pwd=1QOWWuFTTjXQcla2T32GQ5MoyAUHhG.1

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- Assignment 5 Feedback from Hans: inputting files should be directly from the zipped folder, NOT the pathway on your own machine. Always write relative paths, NOT absolute paths
- Research question: Is the microbiome a mediator of the increased prevalence of asthma in patients with MS compared to the general public?
    - More specific question: Determine whether the microbiota differs in RRMS and PMS patients with and without asthma.
- Metadata discussion
    - Variables retained: sample-id, site_x, sex, age, year_of_onset, disease, disease_course, disease_duration, diet_no_special_needs,         diet_special_needs, birth_method, breastfeeding, allergies, allergy_specific, asthma, eczema, ms_family, nsaids, nsaids_specifics,        rxmeds, rxmeds_number, otc_meds, otc_number, smoke
    - Are there any variables missing that would be valuable for our analysis?
          - We haven't included MS treatment variables
    - Are there any variables retained that are not important for the analysis?
          - We wanted to keep variables that would have potential impacts on the microbiome.
          - Hans: Too much, narrow down research question (literature review) to learn about which variables we want to investigate and start writing background and introduction
  - 3 Aims (what we want to do to achieve our research question): Have meaningful link and story from one aim to the next. E.g: Taxnomic analysis of microbiome of MS people vs Non-MS people 
    - Concerns: how do we control for confounding variables?
      2. Determing if there is a difference between those with Asthma and non-asthma
        - Hans: removing confounding variables that do not meet criteria. Run correct statistical analysis. No need to do all variables, only do the ones interested in. Have a targeted approach, choosing a variable from a literature review.
        - Minimum smoking sample size: 2 is too little, more than 10 is ideal
    - There aren't duplicates for every sample - should we filter out the duplicates?
        - Hans: filter out duplicates. 
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
      - Hans: try making manifest file, can refer to ChatGPT with caution
- Hans: when subsetting a sample, set the seed to the same number-> same result (reproducible). Setting one seed should be representative of your whole sample. Subset it once is enough. If you want to check, subset a few times.
- Filter the metadata in R once variables are confirmed
- To Do for next week:
      - Research into the literature. Figuring out objectives and aims, Background and     Introduction
      - Try to make a manifest file
      - Proposed Approach, Overview Flowchart
      - Allocating tasks to each member

Meeting adjourned October 9, 2024 3:36 PM


