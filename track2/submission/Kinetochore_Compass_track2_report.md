# Helping cells cope with chromosome imbalance

## A drug-repurposing proposal for mosaic variegated aneuploidy

**Kinetochore Compass · Track 2 · 17 October 2026**

### Our proposal

A genetic diagnosis tells us where a disease begins. It does not necessarily tell us which consequences we can change. In mosaic variegated aneuploidy (MVA), errors in cell division can leave cells with extra or missing chromosomes. The resulting imbalance may place a continuing burden on the cell's ability to make, fold and remove proteins.

We propose to test whether **sodium 4-phenylbutyrate (4-PBA)** can relieve that burden and improve cell function. We would compare it with **N-acetylcysteine (NAC)** to test the contribution of oxidative damage, and add **felodipine** if impaired protein clearance is demonstrated. These are separate experimental treatments, not a proposed combination.

Our lead is supported by experiments in human cells with chromosome imbalances, including neurons. The important unanswered question is whether those results transfer to MVA at a realistic drug exposure. We therefore define success as improved function after a deficit has appeared, without worsening chromosome errors or favoring the growth of abnormal cell populations. This submission presents a literature-supported hypothesis and a validation plan; drug response has not yet been measured in this case.

### 1. From the genetic finding to a testable mechanism

The leading hypothesis carried forward from our Track 1 analysis identifies **BUB1B**, which encodes **BUBR1**. The table links this proposal to the nominated variant(s) and their recorded sequence annotations. Coordinates use GRCh38 and chromosome:position:reference:alternate notation.

| Variant | Selected transcript | Recorded sequence annotation |
|---|---|---|
| `chr15:40209701:T:G` | ENST00000287598.11 | premature stop signal |
| `chr15:40220612:T:G` | ENST00000287598.11 | amino-acid substitution |

For a two-variant explanation, the saved evidence does not establish that the variants are on opposite chromosome copies. Sequence annotations predict an effect; they do not measure residual protein activity.

BUB1B makes BUBR1, a protein that helps control when chromosomes separate and how they attach to the spindle. Published MVA experiments show that different disease-associated changes can reduce the amount of BUBR1 or leave a protein that does not work properly. Our working interpretation is loss or impairment of this function; the sequence finding alone does not show how much useful protein remains. [Primary mechanism study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2887387/).

**What we need to establish in the cells:** Measure BUBR1 abundance and location, the checkpoint response and errors during completed cell divisions. These measurements distinguish an unstable but potentially useful protein from one whose function cannot be restored simply by increasing its amount.

Our repurposing hypothesis addresses a possible consequence of this defect: a burden of chromosome imbalance that the cells may be able to manage better. We have not measured protein aggregation, oxidative injury or impaired clearance in this case. Those are entry criteria for the experiments below. Any claim that 4-PBA restores the altered disease protein itself would require an additional functional test.

### 2. Why chromosome imbalance could create a treatable problem

Before a cell divides, it copies its chromosomes and distributes them between two daughter cells. The spindle moves the chromosomes; the spindle assembly checkpoint helps delay separation when their attachments are not ready. Defects in this machinery can produce cells with different abnormal chromosome counts.

Two terms matter here. **Aneuploidy** means an abnormal chromosome count. **Chromosomal instability** means that new chromosome-segregation errors keep occurring. A drug might help a cell tolerate an abnormal chromosome count without fixing the process that produced it. Our experiments must measure these outcomes separately.

Chromosomes contain instructions for making proteins. Many proteins work together in assemblies that need balanced amounts of their components. An extra chromosome can increase production of some components without increasing their partners. Missing chromosomes can create shortages. The cell must manage that imbalance while continuing its normal work.

This is where **proteostasis**, the cell's system for maintaining useful proteins, becomes relevant. It includes folding proteins into working shapes, moving them to the right place and removing damaged or surplus material. Human trisomy studies link chromosome imbalance to protein aggregation and cellular stress. They provide a reason to investigate this mechanism in MVA, but do not establish that it is present in every affected cell. [Nawa et al., 2019](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0219592), [Hirata et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7441064/).

Our proposed chain of reasoning is:

**Disrupted cell division → chromosome imbalance → potentially harmful protein stress → a possible opportunity to preserve cell function.**

The final two steps are the hypothesis to test. In particular, Down-syndrome models usually carry a specific extra chromosome, whereas MVA can involve a changing mixture of gains and losses. We are transferring a possible cellular mechanism, not assuming the diseases are equivalent.

### 3. Why we would start with 4-PBA

4-PBA has two advantages for this proposal: the actual compound has been tested in human aneuploid cells, and sodium phenylbutyrate is already an approved medicine for urea-cycle disorders. That existing use provides formulation and pharmacology experience. Its approval does not establish benefit in MVA. [FDA product information](https://www.accessdata.fda.gov/drugsatfda_docs/label/2026/216513s007lbl.pdf).

The most relevant evidence comes from three related studies:

| Study | Finding that supports the proposal | Important limitation |
|---|---|---|
| Human trisomy fibroblasts, 2019 | Reduced protein aggregation; improvement in a senescence-associated measure in a trisomy-21 comparison | A non-neuronal model using millimolar exposure |
| Human trisomy-derived neurons, 2020 | Reduced aggregation, stress and cell-death measures | Treatment began during neuronal differentiation |
| Human trisomy neuronal and glial models, 2026 | Further evidence of aggregation prevention and reduced extracellular phospho-tau | Overlapping investigators; some findings involved a separate genetic intervention |

Sources: [Nawa et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0219592), [Hirata et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7441064/), [Nambara et al.](https://www.nature.com/articles/s42003-026-09902-6).

The 2026 neuronal experiments used 100 µM 4-PBA from the start of neuronal induction through day 21. This supports prevention of certain cellular abnormalities during development. It does not show recovery after an established deficit. A separate co-culture experiment combined 4-PBA with correction of astrocytic DYRK1A gene dosage; that combined result cannot be attributed to the drug alone. The studies also do not constitute independent replication across laboratories. [Nambara et al., methods and Figures 2 and 8](https://www.nature.com/articles/s42003-026-09902-6).

These limitations shape our experiment. We would test treatment after a reproducible deficit appears, and we would require functional improvement rather than a change in aggregate staining alone.

4-PBA is often described as a chemical chaperone—a compound that may help cells handle proteins. Its effects are broader than that description suggests: experimental work also shows altered protein trafficking. We therefore would not assume that fewer aggregates mean the original disease protein has been repaired. [SEC24 trafficking study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5464768/).

### 4. The exposure question must be answered early

A medicine can work in a dish yet fail to reach the relevant cells at a sufficient concentration in the body. This is a central uncertainty for 4-PBA.

The often-cited evidence of entry into cerebrospinal fluid comes from intravenous administration in three nonhuman primates. It does not establish the exposure achieved by an oral formulation in a child's brain. Protein binding, metabolism and exposure duration also matter: a brief peak in total blood concentration is different from sustained drug availability to neurons. [Berg et al., 2001](https://pubmed.ncbi.nlm.nih.gov/11391852/), [TGA assessment](https://www.tga.gov.au/sites/default/files/auspar-sodium-phenylbutyrate-180123.pdf).

We would therefore compare concentration ranges and short versus sustained exposure, with measurement of free drug or intracellular exposure where feasible. If benefit requires an implausible exposure, 4-PBA should lose priority. Laboratory concentrations would not be converted directly into a pediatric dose.

### 5. What the comparison compounds would tell us

**NAC tests whether oxidative damage is a major contributor.** Oxidative stress occurs when reactive molecules overwhelm the cell's defenses and damage its components. NAC can support glutathione, an important part of those defenses. Approved acetylcysteine formulations include an antidote for acetaminophen poisoning. In human Down-syndrome neural cultures, NAC protected against experimentally induced oxidative-stress apoptosis; another trisomy-neuron study found reduced amyloid-beta secretion. These are useful actual-compound results, although neither demonstrates MVA rescue. [Official acetylcysteine label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=472f158a-5ab9-4308-8e49-1116e6ea3d39), [Briggs et al., 2013](https://doi.org/10.1002/stem.1297), [isogenic-neuron study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8405674/).

The evidence is mixed across endpoints. NAC did not reproduce 4-PBA's senescence result in the 2019 fibroblast experiment. Its oral brain target engagement is also unresolved, and a 2026 pediatric NF1 trial did not improve its primary motor outcome. Those differences make a direct comparison more informative than assuming that either protein stress or oxidative stress is the dominant problem. Evidence for N-acetylcysteine amide, or NACA, is not counted as evidence for NAC; they are different compounds. [Fibroblast comparison](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0219592), [oral NAC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5762253/), [NF1 trial](https://pubmed.ncbi.nlm.nih.gov/42173046/).

**Felodipine tests whether improving protein disposal helps.** Autophagy is one route by which cells collect and break down material. Felodipine, an approved antihypertensive medicine, improved selected protein-clearance and disease-related outcomes in neurodegeneration models. Its mouse experiments considered plasma exposure relative to human use, although matching a human plasma peak does not establish equivalent brain exposure. [Siddiqi et al., 2019](https://www.nature.com/articles/s41467-019-09494-2), [official label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=2e4298cf-02e9-49c7-9803-887161e1989a).

We would add felodipine only if the model shows impaired clearance and retains the capacity to improve it. More autophagic structures can indicate either increased disposal or a blocked pathway; a dynamic clearance measurement is essential. Its effects on blood pressure and neuronal calcium signaling also require attention. It has no demonstrated MVA benefit.

| Compound | Question it helps answer | Place in the study |
|---|---|---|
| 4-PBA | Can reducing protein stress preserve function? | Initial comparison |
| NAC | Can reducing oxidative injury preserve function? | Initial comparison |
| Felodipine | Can improving protein disposal preserve function? | Add if clearance is impaired |

We would select the lead from the results. NAC or felodipine could replace 4-PBA if it provides more convincing functional benefit at a feasible exposure.

### 6. Why we did not choose a more direct mitotic drug

Directly improving chromosome segregation is attractive, but the available results require careful interpretation. Very low paclitaxel concentrations reduced ongoing chromosome instability in selected cancer cells with excessive microtubule assembly. However, those cells also grew better in soft agar. This illustrates a central concern: improving an abnormal clone's fitness is not necessarily a desirable outcome. [Ertych et al., 2014](https://www.nature.com/articles/ncb2994).

Other options also depend on a matching mechanism. Rapamycin partially rescued a neural-progenitor phenotype in a fly CIN model, but separate aneuploidy research found that relieving already-depressed protein synthesis could be beneficial. Uridine is another plausible pathway candidate, yet a key chromosome-instability experiment used a mixture of nucleosides rather than uridine alone. These findings support conditional experiments, not a general treatment recommendation. [Fly CIN study](https://www.nature.com/articles/s41467-026-70521-0), [translation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7299149/), [nucleoside-mixture study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4636055/).

Our shortlist favors a small comparison that can answer distinct biological questions. A broader list would not compensate for the absence of a measured, relevant defect.

### 7. The experiment that would decide whether to advance

**First, confirm that the model represents the problem.** We would establish the relevant division defect and measure protein stress, oxidative damage and clearance before selecting drug arms. If stem-cell-derived models are used, we would examine the original cells, multiple independently derived clones and the final differentiated cells. Reprogramming can change chromosome composition or conceal instability; a normal chromosome count at one time does not prove faithful division. [Human reprogramming study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8912248/), [instability study](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002913).

**Second, compare compounds after a deficit appears.** The initial study would include vehicle, 4-PBA and NAC, with appropriate healthy or gene-corrected controls where available. Felodipine would enter the clearance branch. Each compound would have its own exposure range. An early-treatment arm would connect the experiment to published work; a later-treatment arm would test whether an existing cellular problem can improve.

**Third, measure useful function and possible harm together.**

| Question | Measurement | Why it matters |
|---|---|---|
| Do the cells work better? | A prespecified functional endpoint appropriate to the affected cell type; for neurons, activity and network function with maturation controls | Survival or a molecular marker alone does not establish recovery |
| Is the intended mechanism changing? | Protein burden and turnover, completed autophagic clearance, or glutathione and oxidative damage, depending on the drug | Distinguishes plausible mechanism from a nonspecific response |
| Are divisions more or less accurate? | Chromosome errors per completed division in dividing cells, alongside division and death rates | Fewer observed errors can result simply from fewer divisions |
| Is the population changing? | Chromosome composition and clonal behavior over time | Detects loss of vulnerable cells or preferential growth of an abnormal clone |
| Is the exposure realistic? | Concentration–response and exposure-duration relationships | Establishes whether a promising cellular effect has a credible route toward further development |

We would randomize treatment allocation, blind image analysis and repeat independent cultures or differentiations. Multiple images of one culture would not be counted as independent biological replicates. Pilot variability would inform sample size; the primary endpoint, meaningful effect and acceptable adverse-effect margin would be defined before the confirmatory comparison. Results would include effect estimates and uncertainty, not only significance tests.

**The decision rule is straightforward:** advance a compound only when functional benefit is reproducible, the proposed mechanism is supported, exposure is plausible and the genome-stability checks are satisfactory. Aggregate reduction without functional improvement would remain a mechanistic observation. Benefit accompanied by worse chromosome errors or abnormal-clone expansion would stop or substantially narrow development. Failure of all candidates would be an informative outcome.

### 8. Contribution and wider use

The contribution is a way to identify useful rescue while avoiding a misleading result. A drug can lower a stress marker by slowing cells, changing protein secretion or eliminating the most affected cells. Our design asks whether cells function better and whether that benefit persists without an adverse shift in the population.

The same approach can be reused across MVA genotypes and other chromosome-instability disorders. The upstream mechanism determines the model and safety checks; measured downstream defects determine which compounds enter the comparison. This permits different cases to reach different answers.

A later computational extension could combine drug-target networks with the direction of drug-induced expression changes, using public approaches such as SHARP. Such analyses would need tissue, dose and cell-state controls, and must distinguish harmful responses from protective adaptation. No network model has been run for this submission. [SHARP primary paper and software links](https://www.nature.com/articles/s43587-026-01161-8).

We prioritize 4-PBA because of the relevant human-cell evidence, while giving competing mechanisms a fair test. The next advance would be a reproducible improvement in an appropriate MVA model at a credible exposure. That result would provide a stronger basis for further investigation than a drug ranking alone.

### Methods and transparency

The following answers follow the organizer’s Track 2 methods form. They are included here so the Markdown report contains the methods description without requiring a separate Excel upload. The literature assessment runs through 15 October 2026; the submission-format revision is dated 17 October 2026.

#### Team name

Kinetochore Compass

#### Approach

We linked the existing Track 1 finding to its published gene function, then distinguished the upstream division defect from possible downstream protein stress, oxidative injury and impaired clearance. The case-specific variant, working mechanism and evidence limits appear in report section 1. We reviewed primary experiments and official drug information, prioritizing the actual approved compound, relevant human models, functional outcomes and plausible exposure. Sodium phenylbutyrate (4-PBA) is the provisional first candidate because of human trisomy-cell evidence. NAC provides a comparison for oxidative injury; felodipine enters only if impaired clearance is demonstrated. Each drug is a separate experimental arm. The proposed study requires functional benefit after a deficit appears, without worse chromosome errors or abnormal-clone expansion, at a credible exposure. No drug-response experiments have been performed.

#### AI-use disclosure

OpenAI Codex, ChatGPT Plus, GPT-6 Astra; model-improvement/data-sharing was disabled before challenge work, as confirmed by the user. AI assisted with public-literature review, coding and writing.

#### Candidate identification

AI-assisted literature search and source-by-source curation, rather than an automatically fitted drug-response model. No exhaustive drug-target database screen or new expression-signature analysis was run. A local script verifies the existing case-evidence bindings and assembles the case section; it does not predict clinical drug efficacy.

#### Review and curation

The review recorded the model, actual intervention, endpoint, treatment timing, concentration, exposure duration and conflicting evidence. We distinguished prevention from reversal, single-drug effects from genetic or combination rescue, and independent replication from studies with overlapping investigators. We did not substitute NACA evidence for NAC or nucleoside-mixture evidence for uridine alone. Candidate decisions and source audits are retained in track2/research_v4. This is a structured narrative review, not a registered systematic review or an independent clinical expert assessment.

#### Public and restricted sources

Drug-repurposing evidence came from public literature, official drug labels, public trial records and original methods/software. The case mechanism additionally uses the organizer-provided, access-controlled challenge data through the existing local Track 1 analysis. Thus the overall case-based work is not described as public-data-only.

#### Public sources

Primary sources include Nawa et al. (2019), <https://doi.org/10.1371/journal.pone.0219592>; Hirata et al. (2020), <https://pmc.ncbi.nlm.nih.gov/articles/PMC7441064/>; Nambara et al. (2026), <https://doi.org/10.1038/s42003-026-09902-6>; Briggs et al. (2013), <https://doi.org/10.1002/stem.1297>; and Siddiqi et al. (2019), <https://doi.org/10.1038/s41467-019-09494-2>. We also reviewed FDA/DailyMed labels, a TGA assessment, primary pharmacokinetic studies and public trial records. Full source links and limitations are in report sections 3–6 and the accompanying research files. SHARP is discussed only as a possible future extension, not a completed analysis.

#### Proprietary sources

No additional proprietary drug/target database or commercial patient cohort was used. The restricted challenge inputs are identified separately above. Case-derived files remain local; the participant will delete them after the hackathon in accordance with its requirements.

#### Variant-mechanism characterization

Report section 1 contains the locally assembled variant identifiers, selected-transcript consequences, phase limitations and gene-specific working mechanism. We distinguish predicted loss or impairment of function from experimentally established loss-of-function or gain-of-function. Gene-level disease literature supports the pathway interpretation; the sequence finding alone does not measure residual protein function. We propose protein/pathway assays and division measurements to test it. Downstream protein stress, oxidative injury and clearance defects remain hypotheses to measure before using them to select drug arms. 4-PBA is not claimed to repair the disease protein.

#### Time and effort

Several iterative literature reviews and source audits, followed by report preparation. Effort was not prospectively time-tracked, so no numerical person-hour or cost estimate is claimed. The literature assessment is current through 15 October 2026; this submission-format and presentation revision is dated 17 October 2026.

#### Method abstract

Mosaic variegated aneuploidy can leave cells with different abnormal chromosome counts. The resulting protein imbalance may create a reversible source of dysfunction even when the original division defect persists. We propose a mechanism-guided comparison of approved drugs, with sodium phenylbutyrate (4-PBA) as the provisional first candidate.

We linked the existing Track 1 finding to published gene function, separating the inferred upstream mechanism from downstream effects that have not been measured in this case. A structured narrative review assessed actual-compound evidence, model relevance, treatment timing, functional endpoints, exposure and contradictory findings. Human trisomy-cell studies support testing 4-PBA, but do not establish efficacy in MVA, reversal of an established deficit or sufficient oral brain exposure. NAC tests the contribution of oxidative injury. Felodipine is conditional on a demonstrated clearance defect. The drugs would be tested separately.

The proposed validation first confirms model fidelity, then compares prevention and treatment after a deficit appears. Advancement requires reproducible functional improvement with supporting mechanism measurements, credible exposure and no adverse change in chromosome segregation or clonal composition. Division rate and selective cell loss are measured to detect misleading rescue.

The strengths are explicit evidence boundaries, a small informative comparison and joint measurement of function and genome stability. Limitations include cross-disease extrapolation, overlapping investigators in key studies, uncertain brain exposure and absence of case-specific drug-response experiments. The reusable output is a testable decision framework, not a calibrated clinical efficacy ranking.
