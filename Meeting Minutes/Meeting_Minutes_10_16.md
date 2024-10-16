# Team 9 Meeting Minutes 
# Meeting Details:
## Date: 
October 16th, 2024

## Location:
Zoom: https://ubc.zoom.us/j/67359231965?pwd=1QOWWuFTTjXQcla2T32GQ5MoyAUHhG.1

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- Research question: Is the microbiome a mediator of increased prevalence of asthma in patients with MS compared to the general public?
- Review proposed research aims
    1. Investigative comparison of variables composition (alpha and beta diversity, what are the significant differences, functional consequences?) of healthy people as control vs asthma
       - Rephrase this aim
       - Expect to uncover how microbiome varies between our control and treatment groups
       - We retained severity score as one of our variables, so we could potentially look at this. 
    2a. Investigative comparison of variables functionality from first aim
       - Make sure that we are exploring an angle that hasn't been assessed before
       - Investigating prevalence of a pro-inflammatory microbiome. Potentially look at microbes associated with inflammation from the                literature and associate that with a specific niche.
       - Look for microbes associated with pro-inflammatory processes, then assess if it's present in our cohort.
       - Think about approach and analysis that we want to run. We need to have a story in mind that prompted our decision to perform this             analysis.
    2a. DeSEQ2 - look at cross-sectional dynamics in healthy patients vs MS/asthma
    3. Explore further variables (socioeconomic status, education, genetics ms_family) as final figure
       - Functional analysis possible and encouraged.
       - Use picrust --> watch module on Canvas for more information
- Research Overview discussion
  - We can't jump to the conclusion that a decreased abundance of certain strains of bacteria play a protective role in MS; however, we         can compare it to trends that we see in the literature
      - Be careful in discussion of relative abundance
  - see what taxa are correlated with differences then compare to datasets
  - Read paper: Reference frames - Morton et al. 2019 
  - Log of ratio between two taxa then compare log data across samples. This gives you a better understanding of any changes as the ratio        removes sampling biases. This is an interesting approach to relative abundancy analyses.
  - Use picrust for functional analysis
        - Limitations of functional analysis: we don't have trasncript information so we are basing our analysis on taxa data
        - We are not observing, we are predicting which may show functional potential, but not absolute truth
            - We are looking at who could be in the community,
        - Prediction is based on who we determine is in the community since we are only looking at a single region (short read) in the 16s             RNA
            - This is less specific and less resolved. We will most likely only resolve to the genera level
            - We won't know specific strain, which will impact our functional analysis.
        - Only address limitations of picrust if we make an aim that addresses the limitation
- Discuss preliminary data analysis on societal factors 
    - p-value of 0.07 between control and MS population for social variables (education and occupation)
        - is it worth keeping as an aim? No it is not
    - no significance between asthma and mode of delivery 
        - could potentially look at microbiome composition? And longitudinal dynamics
        - possible new aim of microbiome chnages in patients based on mode of delivery and prevalence of asthma and MS
    - very significant association between sexes 
- Metadata Update 
    - right now our metadata has 22 variables (including sample ID) and 152 samples 
- Manifest file progress: able to use the ms_manifest.tsv file from project_2 directory on class server. Is this okay?
- Help with Proposed approach 
- To Do for next week:
      
Meeting adjourned 
