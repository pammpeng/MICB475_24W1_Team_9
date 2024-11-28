# Team 9 Meeting Minutes 
# Meeting Details:
# Date: 
November 27, 2024

# Participants:
Tia Murdoch, Asmita Jain, Annie Saint, Pamela Peng, Claire Rollins

# Agenda Items

## Alpha Beta Diversity (Figure 1)
- Alpha: All metrics except simpsons-> Asthma lowest diversity, MS, MS + asthma, healthy has similar higheset diversity
- Beta: only significant differences in healthy vs MS (p=0.006)
- significant differences in evenness

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

- have to use p-adjusted value not p-value
- not a lot, you can comment on individual genera, as long as you don't ignore the general pattern
- there does not seem to be a big difference

## Picrust2 (Table 1)
- MS_Asthma is most clustered, does that mean they are more similar? less diverse? - just means samples are 

<img width="425" alt="Screenshot 2024-11-25 at 5 36 55 PM" src="https://github.com/user-attachments/assets/b4fb01fc-3aaa-4bc7-83ee-684a62aed455">
- put into supplimental 

- write control_asthma as just asthma
- asthma seems to be drives lots of functional change - matches to alpha diveristy (potential decrease in diveristy, causing decrease in functionality)
- go through all the pairwise comparisons and tabulate the results - gan gauge which comparisons to keep the plot of
- p =0.05, log +/- 1.5
- if you reference specific pathways then add it to the supplimental
- S1 - PCoA graph 

## Proposed Narrative
- MS is potentially masking the effects of asthma due to similarities in MS+Asthma vs MS comparisons
https://docs.google.com/document/d/1yV2v-xm3Y2zUjSDU6DhI1xWOvH8VZ7FDswyeJf18VB8/edit?usp=sharing

- not all lot different at a taxa level, but more differences in functional analysis

Meeting Monday 9:30 am for feedback on slides, keep Wednesday meeting 
