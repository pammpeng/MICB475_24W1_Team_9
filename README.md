# MICB475_24W1_Team_9_Meeting_Minutes

## Date: 
October 2, 2024
## Participants: 
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng
## Absentees: Claire Rollins

## Agenda Items
1. Research Area Ideas:
  - Nasal Microbiome and Asthma and Smoking in different countries, socioeconomic status, living conditions
   - MS dataset with Asthma
   - 2 Control Groups: 40 without Asthma and MS, 40 with MS and without Asthma.
   - Treatment group: 40 with Asthma and MS.
   - Technical replicates (biologically different).
   - Three types of MS: RRMS correlates with Asthma only, very little Primary MS and Secondary MS don't have Asthma, just do RRMS
   - Use R to wrangle data on local computer (manifest) then make new files for server (qiime)
   - Summarize sample sizes first (tabulate how many sample sizes for each condition)
   - Put all code on Github now, annotate code properly with comments!!
3. Brainstorm/narrow down research question
  - Is the microbiome a mediator of the increased prevalence of Asthma in patients with MS compared to the general public?
4. Meeting online next week!

## Date: 
October 9th, 2024

## Location:
Zoom

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

## Date: 
October 16th, 2024

## Location:
Zoom

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

## Date: 
October 23rd, 2024

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- Research question: Does the gut microbiome mediate the increased prevalence of asthma in patients with multiple sclerosis (MS)?
- Research aims:
    1. Investigate Alpha and Beta Diversity in the Gut Microbiome of MS patients with and without Asthma
    2. Identify Differentially Abundant Microbial Taxa Using DESeq2 Analysis and Core Microbiome Identification
          analysis. Hans: Log Ratio is a better approach; Ancom and Aldex2 tools, plug ins into QIIME2
    3. Conduct Functional Analysis to Explore Microbial Metabolic Pathways using PICRUSt2
       - Module 17: PICRUSt2 as a resource, use internet, look at previous MICB475 GitHub codes just reference it
- Proposal debrief
    - What went well?
    - What could be improved for next time?
- Tasks already completed
    1. QIIME2 Processing
- To-do for October 29th
    1. Alpha and Beta diversity analysis
    2. Taxonomic analysis-> requires literature search
    - Start this analysis this week ideally so it is ready.
    - Tip: Start writing results and methods as we go 
- Assign roles
## Date: 
October 30, 2024

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng

## Agenda Items
- Proposal debrief and revision based on comments, receive feedback from Hans
- Submit version November 5th Tuesday. 
- Should we revise our aims?
- Tasks already completed
    1. QIIME2 Processing
- To-do
    1. Alpha and Beta diversity analysis-> no significant differences in alpha diversity p value= 0.189, significant difference for beta diversity (unweighted unifrac) in sex, try pCoA.
    2. Taxonomic analysis-> requires literature search
Assign roles
- DESEQ2 -> Annie and Tia. Hans suggests AnCom or Aldex2 (log ratios) over DESEq2.
- PICRUSt2 -> Asmita and Pam

- Report all code with comments on each code.
  ## Date: 
November 6, 2024

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- Done by Tia
  1. Alpha and Beta diversity 
    
- To-Do
  1. DESEQ2 -> Annie and Tia. Hans suggests AnCom or Aldex2 (log ratios) over DESEq2.
  2. PICRUSt2 -> Asmita and Pam Updates. Issues with R visualizations

Notes:
- Presentations on Dec 3 and 5
- Only 2-4 Meetings left!
- Literature review for taxonomy
- Optional Meeting next week. Wednesday 11-2 or after 2:45pm, Thursday before 10:15 or after 12 till 2pm.
- Put codes on GitHub
- Use ms_metadata_final_2
  ## Date: 
November 20, 2024

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- PiCrust2 interpretation feedback
- p<0.005, and filter out fold change +/- 1.5
- Find pathways of interest in literature 
- DESEq2 analysis feedback
-   Ancom help, code failing (error= null)
-   Watch Ancom tutorial
-   make summary table with up and downregulated genes, add picrust 2 deseq
- Alpha beta diversity Results by Tia
-   krusko wallis okay
-   Asthma only affects microbiome if you don't have MS
-   MS already depletes microbiome enough
-   Make axis labels bigger
-   Don't need group label
-   Having asthma doesn't affect MS
-   show variants
-   Add ellipses anyways
-   also do for unifrac, weighted unifrac
-   Axis lables need to be bigger
-   Colours for taxonomy plot are not great
-   remove ticks (ask ChatGPT)
-   Start thinking about how you will format figures. Add * and p values.
- Core Microbiome by Claire

Conclusions: 
- Asthma affects only healthy patients
- MS masks effects of all other aspects like asthma
- Make summary key points for each figures

# Date: 
November 27, 2024

# Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

# Agenda Items

## Alpha Beta Diversity (Figure 1)
- Alpha: All metrics except simpsons-> Asthma lowest diversity, MS, MS + asthma, healthy has similar higheset diversity
- Beta: only significant differences in healthy vs MS (p=0.006)
- significant differences in evenness

   **Notes**
- remove title from panel A
- use ggsignif for signifiance on graphs
- asmtha decresease alpha diversity
- don't think we need to show c
- don't bother trying to format nicely
- a: remove title, make axis labels bigger, numbers are too small, want to include breaks that are above and below highest and lowest points (y-axis), include 1 and 5, want to scale the figures the same size
- b: change the theme to be the same as a (think a is classique), would be nice in the box plot if the colours were the same as pannel b
- could include the results from c to the figure ledgend, no signifiance marker on PCoA,
- figure ledgend: move away from talking about it in the technical sense, talk about it in a biological level, add significance threshold
- Takeaways: asthma effects the alpha diversity but MS effects the beta diversity, beta diversity is more sophisticated, don't know what to make of these results (Evelyn), doesn't trust the p-value for beta diversity, no visible clustering, trust panel B 

## Core microbiome (Figure 2)
- control group has the most unique ASVs
- highest degree of overlap in composition between MS and MS_Asthma
- control vs MS has the same overlap as control vs asthma
- MS_Asthma vs Control = 62%
- Asthma vs Control = 71%
- MS vs Control = 71%
- MS vs MS_asthma = 78%
- **overall MS masks the effects of asthma but asthma still plays a role?**

**Notes**
- share more of a core, pattern shows there is no significance between the groups, don't think the individual venn diagrams add
- core microbiome is mostly shared 

## DESEQ (Figure 3)
- MS vs Control: 15 decrease, 8 increase
- Asthma vs Control: 27 decrease, 6 increase
- MS+Asthma vs Control: 19 decrease, 9 increase
- MS+Asthma vs MS: 13 decrease, 10 increase
- MS+Asthma vs Asthma: 5 decrease, 18 increase
- MS vs Asthma: 8 decrease, 19 increase
- Perform lit search to investigate role of specific genera?

**Notes**
- look at more patterns, reference on right, tabulate
- asthma majority decreases matches alpha diversity
- decrease is generally the pattern, there is something, both conditions generally decrease the taxa but compared to each other, MS is increasing relative to Asthma
- 4 relevant pannels (vertical) - MS vs Control, Asthma vs Control, MS+Asthma vs Asthma, MS vs Asthma
- go through the bars and colour code all the genera that are consistant (more of a rational to pick though genra to talk about)
- non-consistant = grey, consistant = orange or red 

## ALDEx2
- Blautia more abundant in asthma vs control group, but less abundant in MS vs Asthma and MS+Asthma vs Asthma
- Romboutsia is even less abundant in MS+Asthma vs Asthma compared to MS vs Control and MS+Asthma vs Control
- No graph generated for MS_asthma_vs_MS as there were no genera that were significantly different between the two groups

**Notes**
- have to use p-adjusted value not p-value
- not a lot, you can comment on individual genera, as long as you don't ignore the general pattern
- there does not seem to be a big difference

## Picrust2 (Table 1)
- MS_Asthma is most clustered, does that mean they are more similar? less diverse? - just means samples are 

<img width="425" alt="Screenshot 2024-11-25 at 5 36 55 PM" src="https://github.com/user-attachments/assets/b4fb01fc-3aaa-4bc7-83ee-684a62aed455">
- put into supplimental 

**Notes**
- write control_asthma as just asthma
- asthma seems to be drives lots of functional change - matches to alpha diveristy (potential decrease in diveristy, causing decrease in functionality)
- go through all the pairwise comparisons and tabulate the results - gan gauge which comparisons to keep the plot of
- p =0.05, log +/- 1.5
- if you reference specific pathways then add it to the supplimental
- S1 - PCoA graph 

## Proposed Narrative
- MS is potentially masking the effects of asthma due to similarities in MS+Asthma vs MS comparisons

- not all lot different at a taxa level, but more differences in functional analysis

Meeting Monday 9:30 am for feedback on slides, keep Wednesday meeting 

## Date: 
December 4, 2024

## Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

## Agenda Items
- feedback on formal figures 
- feedback on final narrative/title
## Notes on Story
- core microbiome and diversity metrics there isn't a difference in diveristy between the two groups
- DeSeq doesn't align with this - MS and asthma alone reduce a lot of genra, add MS to asthma restores some of the genera 
- Picrust didn't get much
- MS is restoring diversity, seen also in the DeSeq
- presenting the DeSeq as a ratio of increase to decrease
- MS seems to have a restoring effect on the decrease caused by asthma 
 - discussion the contriditing evidence in the discussion section
 - we are leading towards the side that there is no influence on the microbiome
 - can still discussion the shared genera
 - still discussion the pie crust, there is also no functional differences 
 - diseases effect the microbiome in different ways
 - could be balancing each other out with functional redundency
 - ## there is not an effect on the microbiome
 - Pam gets to include her PcoA plot!!!!!! (cause it shows there is nothing)
 - move the table to sup and put figure 4 as PcoA plot

