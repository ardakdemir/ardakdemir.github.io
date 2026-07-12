---
layout: post
title: "My Takeaways from ACL 2026"
description: "Notes and reflections on the research themes that stood out to me at ACL 2026."
---
# Working Drafts

This document contains both versions while they are being developed. Remove the LinkedIn section before publishing the blog post.

## Version 1: LinkedIn Post

I just got back from ACL 2026 and thought I can share a few things that stuck with me.

My main takeaways:

- **The submission crisis is a “success catastrophe.”** One point I took from Philip Resnik's keynote was that NLP has grown so much that reviewing can no longer keep up. I also came away thinking about how similar ACL, EMNLP, and even broader ML conferences have become, and the growing tension between publishing to advance knowledge and publishing to compete for top research jobs.
- **Evaluation does not just measure progress; it shapes it.** Barbara Plank made a strong case for taking evaluation more seriously as a science and for keeping language at the center of NLP.
- **Build work that compounds.** Comments in the explainability panel and the Big Picture invited talk pushed back on the short-term publishing mindset: jumping from one small, publishable gain to the next instead of doing deeper, long-lasting research.
- **Strong researchers are also strong engineers.** In the Careers in NLP panel, Iz talked about how blurry the line between research and engineering has become. Taking ownership, executing, and being able to unblock yourself are all important research skills.
- **Poster halls beat oral sessions—for me.** Posters made it much easier to find unexpected work, meet people, ask questions, and jump between topics. I got a lot more out of that format.
- **NLP is still incredibly broad.** Even though so much of the field is converging around LLMs, I saw a huge range of ideas and research questions—especially some really interesting work grounded in language and linguistics.

Despite all the challenges the field is facing, I came back super energized by how many interesting ideas there are still to work on. Exciting times for NLP! And finally, San Diego weather was amazing :)

I wrote a longer version of my takeaways here:
https://ardakdemir.github.io/blog/2026/07/12/acl-2026-takeaways

---

## Version 2: Detailed Blog Post

## Introduction

I just got back from ACL 2026 and wanted to write down the talks, papers, and conversations that stayed with me. I left with some concerns about where NLP is heading: too many submissions, conferences becoming harder to tell apart, short-term publishing incentives, and more of the field converging around LLMs.

At the same time, I came back very inspired. I met researchers working on all kinds of problems, heard people argue for deeper and more durable research, and saw a lot of creative work. These are not complete conference notes—just the things I found most interesting.

## Overall Insights

- **NLP is in a submission crisis.** Philip Resnik called it a “success catastrophe”: the field has grown faster than the reviewing system can handle, and no one seems to have a full solution. The review crisis is really a symptom of this bigger problem. See this [discussion of the ACL 2026 business meeting](https://www.linkedin.com/posts/nikos-aletras-6b797422_the-acl-2026-business-meeting-during-the-activity-7480159920043102208-eLHb?utm_source=share&utm_medium=member_desktop&rcm=ACoAABOiqN8BLg2C5pbHHEjJHpMXDOJLmZ7F1T4).
- **Build scaffolds, not just papers.** Codebases, datasets, workshops, and communities give future work something to build on. A single paper can produce a result; a community can build a field. Noah Smith, Tal Linzen, and others came back to this idea.
- **Evaluation shapes what we build.** Barbara Plank argued that evaluation should be treated as a measurement science, not as a final box to check. What we choose to measure ends up shaping what the field works on.
- **Keep language central to NLP.** Many of the papers I found most interesting were grounded in language or linguistics instead of treating every problem as another generic LLM application.
- **Research and engineering are getting harder to separate.** In the Careers in NLP panel, Iz argued that good researchers also need to be good engineers. They take ownership, unblock themselves, and get things done.
- **Poster sessions were a much better use of my time.** There was more interaction, more flexibility, and more chance to stumble onto unexpected work. The posters were at least as interesting as the oral papers, and it was easy to move around and meet people. Several oral sessions I attended also lost time to small issues with connections or slides.
- **The field is broader than it first looks.** Even with so much attention going to LLMs, I was pleasantly surprised by the variety and creativity of the research. [Your Students Don't Use LLMs Like You Wish They Did](https://aclanthology.org/2026.acl-long.875.pdf) was one memorable example.
- **Brand still matters.** Posters and authors from top institutions clearly attracted more attention. The crowd size did not always reflect the quality of the work.

## Highlights from Keynotes and Talks

### Panel Discussion on Explainability

It was brought up during the panel that academia rewards small, publishable findings even when they do not add up or generalize. For explainability research in particular, we need deeper work that holds up across different models, domains, and use cases.

This came up again and again during the conference: do not keep jumping from one small improvement to the next. Pick questions and build work that can last.

### Barbara Plank's Presidential Address

Two points from Barbara Plank's address really stuck with me. First, language should stay central to NLP research. Second, evaluation does not just report progress—it shapes progress by deciding what the community optimizes for.

Her suggestion was to put much more focus on **evaluation as a measurement science**. This feels especially important for LLM-as-a-judge methods, where there is still a lot we do not understand.

### Opening Keynote by Philip Resnik

Philip Resnik described the submission crisis as a “success catastrophe.” NLP's growth is a success, but the scale has overwhelmed the systems the field depends on. The review crisis is not the root problem; it is the most visible symptom.

He also pointed out how hard it has become to tell ACL and EMNLP apart. ACL used to represent computational linguistics more broadly, while EMNLP had a clearer focus on empirical and statistical methods. Today, both are dominated by similar transformer-based work, and even NeurIPS and ICML can feel increasingly similar.

That is worrying because science needs people to try different paths, including paths that fail. One possible response he mentioned was taking venue fit and the type of contribution more seriously during review.

His most uncomfortable point was about the imbalance between publishing to advance human knowledge and publishing to get jobs at top industry labs. Right now, the incentives lean heavily toward the latter.

### Noah Smith at the Big Picture Workshop

In his invited talk, [“Where It Hurts: Finding Durable Questions While Moving Fast,”](https://www.bigpictureworkshop.com/) Noah Smith used a trigger-point analogy for choosing research problems. The idea is not to chase every immediate pain point, but to find the places where focused work could unlock broader progress.

He contrasted two failure modes:

1. **Symptom chasing:** working on timely, short-term questions that produce papers but never add up to a lasting research program.
2. **Protective stiffness:** ignoring the rest of the field and focusing only on one durable question until the work becomes isolated in a very small niche.

The useful middle ground is to find high-leverage pressure points: maybe an evaluation method everyone knows is limited but keeps using, or a dataset the field needs but no one has built. Smith's [publication history](https://nasmith.github.io/publications/) gives more context for the kind of durable research he had in mind.

### Careers in NLP Panel

Iz had a very practical take on the blurry line between research and engineering. Successful researchers do not stop at ideas. They take ownership of the whole effort, build what they need, unblock themselves, and execute.

I really liked this framing. In modern NLP, engineering is not just support work around research; it is often part of the research itself.

## Selected Papers

### AutoRAN: Automated Hijacking of Safety Reasoning in Large Reasoning Models

- **Authors / venue:** Jiacheng Liang, Tanqiu Jiang, Yuhui Wang, Rongyi Zhu, Fenglong Ma, and Ting Wang (Stony Brook University and Penn State) — ACL 2026 Long Paper, board G189.
- **What it is about:** AutoRAN is the first framework to automate jailbreaking of large reasoning models by attacking their safety reasoning traces. A weaker, less-aligned model simulates the target's execution reasoning, wraps a harmful goal in an educational narrative, and iteratively refines prompts using reasoning leaked through the target's refusals. It reaches nearly 100% attack success within one or a few turns on GPT-o3/o4-mini and Gemini 2.5 Flash across AdvBench, HarmBench, and StrongReject. The key takeaway is that reasoning transparency can itself become an attack surface.
- **Links:** [arXiv](https://arxiv.org/abs/2505.10846) · [ACL Anthology](https://aclanthology.org/2026.acl-long.1988/) · [Code](https://github.com/JACKPURCELL/AutoRAN-public)

### The Agent's First Day: Benchmarking Learning, Exploration, and Scheduling in Workplace Scenarios

- **Authors / venue:** Daocheng Fu, Jianbiao Mei, Rong Wu, Xuemeng Yang, Jia Xu, Ding Wang, Pinlong Cai, Yong Liu, Licheng Wen, and Botian Shi — ACL 2026 Findings, boards G30/G186 area.
- **What it is about:** TraineeBench is a dynamic agent benchmark that simulates a trainee's first day at work. It evaluates three capabilities missed by static benchmarks: context-aware scheduling of streaming, prioritized tasks; active exploration under uncertainty; and continual learning from previous tasks. Frontier agents show substantial gaps, particularly in exploration and continual learning.
- **Links:** [arXiv](https://arxiv.org/abs/2601.08173) · [ACL Anthology](https://aclanthology.org/2026.findings-acl.1505/) · [Code (EvoEnv)](https://github.com/KnowledgeXLab/EvoEnv)

### Speculative Verification: Exploiting Information Gain to Refine Speculative Decoding

- **Authors / venue:** Sungkyun Kim, Jaemin Kim, Dogyung Yoon, Jiho Shin, Junyeol Lee, and Jiwon Seo (Hanyang University and Seoul National University) — ACL 2026 Findings, poster #5555-FIND.
- **What it is about:** Speculative decoding wastes compute when the draft model is wrong, especially with large batches. Speculative Verification adds a small companion model that estimates draft–target alignment through information gain and dynamically adapts verification length. It requires no changes to the draft or target models and reports up to 1.9× speedup, averaging roughly 1.4× at large batch sizes, across 13B–72B target models and seven NLP tasks.
- **Links:** [arXiv](https://arxiv.org/abs/2509.24328) · [OpenReview](https://openreview.net/forum?id=k2DaAIB55l)

### DQA: Diagnostic Question Answering for Enterprise IT Support

- **Authors / venue:** Vishaal Kapoor, Mariam Dundua, Evren Yortucboylu, Sarthak Ahuja, Neda Kordjazi, Yiming Li, Vaibhavi Padala, Derek Ho, Jennifer Whitted, and Rebecca Steinert (Amazon) — ACL 2026 Industry Track, board C25.
- **What it is about:** DQA treats IT support as a diagnostic problem in which the root cause is latent, rather than as plain question answering. It maintains a persistent diagnostic state and uses Retrieval-Aggregated Generation to aggregate retrieved cases at the root-cause-cluster level. Across 150 enterprise scenarios, it reports a 78.7% success rate versus 41.3% for multi-turn RAG, while reducing the average number of turns from 8.4 to 3.9.
- **Links:** [arXiv](https://arxiv.org/abs/2604.05350) · [ACL Anthology](https://aclanthology.org/2026.acl-industry.79/) · [Amazon Science](https://www.amazon.science/publications/dqa-diagnostic-question-answering-for-it-support)

### Toward Scalable Verifiable Reward: Proxy State-Based Evaluation for Multi-turn Tool-Calling LLM Agents

- **Authors / venue:** Yun-Shiuan Chuang, Chaitanya Kulkarni, Alec M. Chiu, Avinash Thangali, Zijie Pan, Shivani Shekhar, Yirou Ge, Yixi Li, Uma Kona, Linsey Pang, and Prakhar Mehrotra (PayPal AI) — ACL 2026 Industry Track, board C29.
- **What it is about:** Evaluating multi-turn, tool-calling agents typically requires an expensive deterministic backend. Proxy State-Based Evaluation instead uses an LLM state tracker to infer structured state from an interaction trace, then uses LLM judges to verify goal completion and detect tool or user hallucinations against scenario constraints. It produces stable, model-differentiating rankings, generates on-policy and off-policy training data, and reports more than 90% agreement between humans and judges.
- **Links:** [arXiv](https://arxiv.org/abs/2602.16246) · [ACL Anthology](https://aclanthology.org/2026.acl-industry.87/)

### Frontier Models Fail at Moderating Pluralistic Social Media

- **Authors / venue:** Zoher Kachwala, Bao Tran Truong, Rasika Muralidharan, Haewoon Kwak, Jisun An, and Filippo Menczer (Indiana University, Observatory on Social Media) — ACL 2026 Long Paper, board C64.
- **What it is about:** PluRule is a multimodal, multilingual benchmark for moderation in pluralistic communities, each with its own rules. Given a comment and its context, a model must identify which rule, if any, was violated. Even GPT-5.2 with high reasoning reaches only about 58%, barely above the 50% “no violation” baseline, while open-weight models fail to beat the baseline. Models handle universal rules such as civility more reliably than context-dependent rules involving effort, evidence, or relevance.
- **Links:** [arXiv](https://arxiv.org/abs/2605.17187) · [ACL Anthology](https://aclanthology.org/2026.acl-long.1590/) · [Code](https://github.com/Zoher15/PluRule)

### Demystifying Multi-Agent Debate: The Role of Confidence and Diversity

- **Authors / venue:** Xiaochen Zhu, Caiqi Zhang, Yizhou Chi, Tom Stafford, Nigel Collier, and Andreas Vlachos (University of Cambridge Language Technology Lab) — ACL 2026 Findings, board C13.
- **What it is about:** Vanilla multi-agent debate often fails to outperform majority vote because homogeneous agents with uniform belief updates preserve expected correctness. The authors propose diversity-aware initialization and confidence-modulated debate. Their theoretical account frames debate as a martingale, while experiments across six reasoning QA benchmarks outperform both vanilla debate and majority vote, especially on harder tasks.
- **Links:** [arXiv](https://arxiv.org/abs/2601.19921) · [ACL Anthology](https://aclanthology.org/2026.findings-acl.1694/) · [Code (DMAD)](https://github.com/SpaceHunterInf/DMAD)

### Evaluation Pitfalls and Sparsity Limitations in LLM-Based Confidence Estimates for Classification

- **Authors / venue:** Elena Merdjanovska, Omar Zaidan, and Andreas Rücklé (Amazon) — ACL 2026 Findings.
- **What it is about:** Verbalized LLM confidence scores are extremely sparse; for example, Qwen3-32B produces only eight unique values on SST-2, with more than half exactly 95%. This sparsity makes evaluation sensitive to how the area under the accuracy-rejection curve is interpolated and can even reverse method rankings. The authors recommend standardized stepwise interpolation and propose verbalization log probabilities, which improve AUARC by 2.3 points without additional inference cost.
- **Links:** [ACL Anthology](https://aclanthology.org/2026.findings-acl.1671/) · [Amazon Science](https://www.amazon.science/publications/evaluation-pitfalls-and-sparsity-limitations-in-llm-based-confidence-estimates-for-classification)

### Statistically Reliable LLM-Based Ranking Evaluation via Prediction-Powered Inference

- **Authors / venue:** Abhishek Divekar and Anirban Majumder (Amazon Science) — ACL 2026 extended abstract; full paper at AAAI 2026.
- **What it is about:** PRECISE uses Prediction-Powered Inference to combine a small human-labeled set with a large LLM-judged set, producing provably unbiased estimates of ranking metrics while correcting LLM-judge bias. It makes hierarchical metrics such as Precision@K tractable by reducing the output space from `O(2^|C|)` to `O(2^K)`. On ESCI, 30 human labels plus Claude 3 Sonnet reduce the Precision@4 standard error by about 21%. In production, 100 human labels were enough to select the best of three variants, later confirmed by an A/B test.
- **Links:** [arXiv](https://arxiv.org/abs/2606.05308) · [AAAI paper](https://doi.org/10.1609/aaai.v40i47.41427) · [Amazon Science](https://www.amazon.science/publications/precise-reducing-the-bias-of-llm-evaluations-using-prediction-powered-ranking-estimation)

### Sycophancy Negatively Affects LLM-as-a-Judge in Conflict Evaluation

- **Authors / venue:** Naghmeh Farzi, Laura Dietz, and Samuel Carton (University of New Hampshire) — GEM Workshop at ACL 2026.
- **What it is about:** The paper tests whether relabeling one speaker as the first-person “Me” biases an LLM judge even when the evidence is unchanged. Across three tasks, four models, and visible and hidden conditions, models systematically favor the narrator and assign less blame or responsibility to the “Me” speaker. The result warns against using LLMs to moderate first-person conversational data without accounting for narrator bias.
- **Links:** [ACL Anthology](https://aclanthology.org/2026.gem-main.45/) · [PDF](https://aclanthology.org/2026.gem-main.45.pdf)

### Other Paper Bookmarked for Later

- **Social Impact Award winner:** [Your Students Don't Use LLMs Like You Wish They Did](https://aclanthology.org/2026.acl-long.875/)

## Conclusion

I left ACL 2026 with mixed feelings. NLP is struggling with scale, incentives, and more of the field clustering around LLMs. But the conference also reminded me how broad and creative the field still is.

The work I want to follow—and hopefully contribute to—is work that compounds: durable questions, language-grounded research, better evaluation, useful infrastructure, and communities that can build more than one paper. I also came away convinced that doing this well takes both research judgment and engineering ownership.

Most of all, I came back energized. Despite all the challenges, there are still so many interesting problems to work on and so many thoughtful researchers approaching them from different directions. Exciting times for NLP!

