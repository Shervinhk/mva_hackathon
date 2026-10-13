# AI use and privacy for this project

Updated 9 October 2026. This document contains workflow instructions and public policy information, not patient findings. It distinguishes the organizer's permission to use qualifying services from the user's stricter instruction to keep patient information local.

## Current boundary

**Patient data and patient-derived results must remain on the local computer.** Do not put them into chat, tool responses visible to an AI model, browser queries, hosted annotation requests, uploaded files, shared conversation links, external repositories, or challenge submission forms. Do not ask an AI assistant to preview a private report or CSV. Local storage alone does not make an AI-assisted workflow local: content returned to a hosted assistant becomes part of its input.

Public literature, public reference resources, source code without patient-specific content, and fabricated test fixtures can be researched and discussed. Deterministic local programs can process the private inputs and save private outputs. Capture their stdout and stderr to local files; return only a generic completion/failure status to the assistant. Review errors locally before sharing a sanitized description. Do not print a patient's allele, phenotype, candidate gene list, score, or clinical narrative as a debugging aid.

On 8 October, macOS reported that this workspace's old `outputs/` folder was included in cloud sync. A new `work/private_analysis.nosync/` folder was created and its exclusion was verified through File Provider metadata before patient execution. The v3 ranker, comparator, wrapper logs and review package now use that folder by default. This check applies to the current computer and folder; check any custom destination separately. Existing inputs and historical outputs require a separate review of their sync and retention state. This change cannot establish whether earlier copies were synced or remove them.

Prepare a local submission package, but do not upload or publish it while this restriction remains in force. A public reference download is permitted; a query constructed from patient findings is not. Avoid implicit context collection too: turn off IDE context when patient files are open, and keep those files away from any active computer-history, screenshot, connector or repository-sharing feature. The official documentation describes how IDE context includes open files. [OpenAI developer settings](https://learn.chatgpt.com/docs/developer-settings).

## What is known and what cannot be guaranteed

Patient-specific information had already appeared in the conversation before the stricter local-only instruction. A new instruction cannot withdraw earlier content from the conversation or establish how it was previously retained or used. Do not claim that patient information has never reached OpenAI, that earlier content was deleted, or that all work was conducted entirely offline.

The source review also encountered case-specific input filenames embedded in earlier script defaults. On 9 October those defaults were replaced with a generic resolver; the actual paths are now in excluded local configuration, and a local check confirmed they are absent from the refreshed source archive. This change does not remove filenames from earlier assistant-visible source or conversation records. Disclosure must distinguish source metadata exposure from locally processed genomic contents and findings; it must not claim that all assistant-visible context after the restriction was free of case-specific information.

The assistant cannot independently guarantee model-training exclusion, zero retention, or the user's account settings. The participant has stated that Improve the model for everyone was disabled long before this work. This is a participant-reported setting, not an independently inspected control. The exact product tier, sign-in route, applicable grant terms, feedback history, and any safety-review status remain unverified. No account settings were changed as part of writing this document. A privacy promise in a prompt is not a substitute for the provider's actual controls and terms.

For future work, a clean chat containing only public code and synthetic examples reduces the risk of carrying earlier private context forward. It does not erase the previous chat. Review account-side deletion and retention options yourself if needed; do not assume clearing a local file removes remote copies.

## Organizer guidance

The organizer text supplied by the user distinguishes a service acting only as a processor from a recipient allowed to use the data for its own purposes. Qualifying services must not train on inputs or outputs or take rights for their own use. Retention must be limited in purpose and time; zero retention is not required. Feedback, safety-review exceptions, and special credit-program terms can change the analysis. The methods report must state the real provider, tier and relevant settings.

The latest public correction in the same discussion permits consumer tiers subject to training opt-out and its stated feedback/safety-review conditions. That permission does not override this project's local-only instruction, and it does not verify this account. [Organizer clarification and subsequent correction](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2).

## Check the actual OpenAI product settings

1. Confirm the active product, account/workspace, tier, and sign-in route in the product interface. Codex can use either ChatGPT subscription sign-in or an API key. The applicable policies differ; the API's default no-training policy must not be assumed to describe a ChatGPT consumer session. [Authentication](https://learn.chatgpt.com/docs/auth), [API data controls](https://developers.openai.com/api/docs/guides/your-data).
2. Review the current consumer data controls in the signed-in ChatGPT account and ensure training opt-out is enabled. The retrieved Codex documentation explicitly links ChatGPT data controls, but does not establish the precise current consumer menu path or the setting's state for this account. Follow the Data Controls FAQ linked from the official page instead of treating an unverified menu instruction as proof. Record the date, setting label and value privately. [Official documentation linking ChatGPT data controls](https://learn.chatgpt.com/docs/customization/computer-history#what-data-gets-shared-with-openai).
3. Do not submit thumbs-up/down ratings or other feedback involving private work. Supported Codex local clients can disable the feedback command with the configuration below. This is an additional feedback control; it does not itself opt out of model training or stop inference traffic. Do not overwrite an existing configuration blindly. [Codex feedback controls](https://learn.chatgpt.com/docs/config-file/config-advanced#feedback-controls).

```toml
[feedback]
enabled = false
```

These checks document actual conditions. They do not authorize patient uploads under the current instruction, and no setting should be recorded as enabled until the participant has verified it.

## Retention and deletion

The organizer's supplied clarification requires deletion of material carrying the genome: raw sequencing and alignment/variant files, indices, copies and subsets, genotype-scale intermediate tables, caches/notebook state, local prompts or logs containing blocks of variant data, and weights or embeddings trained directly on the raw genomic data. Include backups and controlled cloud locations in the inventory. The general rules specify deletion within 30 days after hackathon close and a deletion confirmation to the organizers; no deletion or email has been performed by this documentation task. [Official rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py).

The clarification permits keeping limited findings such as the ranked submission, HPO terms, gene/pathway analyses, and code/reports without prohibited genomic material. A small named finding differs from a table that reconstructs a meaningful part of the genome. Provider logs outside the participant's control are addressed through the provider's permitted terms; the participant should attest only to systems they control. Keeping findings does not mean this user has authorized publishing them. [Organizer clarification](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2).

## Methods disclosure to complete locally

The following is a draft, not a completed attestation. Replace every placeholder with verified facts before submitting any report. If a fact cannot be established, say so rather than inserting a favorable assumption.

> OpenAI Codex was used for code development, public-source research and documentation. Product and sign-in route: [TO CONFIRM]. Product plan/tier: [TO CONFIRM]. Model(s) and dates of use: [TO CONFIRM]. Training-control label, value and date verified: [TO CONFIRM]. Feedback history and relevant credit-program terms: [TO CONFIRM]. Patient-specific content had entered the conversation before the local-only restriction was set. Later source review also encountered case-specific input filenames; these were subsequently removed from the shareable code. The v3 genomic processing and findings were handled by local programs with detailed output kept in excluded local files. No retrospective no-training or no-disclosure guarantee is claimed. Actual implementation and compliance of these boundaries: [TO CONFIRM BY PARTICIPANT REVIEW].

Describe which local programs produced the submission and which parts received human review. Disclose earlier exposure to public participant methods where relevant to scientific independence. The submission validator checks formatting and unresolved placeholders; it cannot certify provider settings, scientific correctness, privacy compliance, or successful submission.
