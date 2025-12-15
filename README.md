FluidBench
================

The goal of FluidBench is to benchmark various LLMs on fluid dynamic
simulations with R. Showcasing the LLMs ability to apply a domain most
likely in the training set, in a programming language not often used for
fluid dynamics.

## FluidBench v1

- bench: five prompts, in-series, fluid simulation tasks
- marks: 1) executable code, 2) solves prompt, correctness, 3) style

### FluidBench Folder Structure

- /prompts contains all questions to ask LLM (01-05) in series

  - the entire .md context is shared with the LLMs one-by-one
  - LLMs are expected to one-shot
    - the LLMs do not receive feedback between prompts

- /bench contains company/model with answers saved in 01-05 files.

  - models return code and/or text with a code snippet - code is saved
    in \*.R while the entire return is saved in .md - if there is only a
    .R then that is all the model returned (e.g. Opus) - code is
    formatted with air for all LLM answers

- /R contains files that *marks* the *bench*. - /bench has the
  standalone answers to the prompts - these are graded with /R

- /data contains results from sourced scripts

### FluidBench Leaderboard v1:

- 1)  can the code be executed without errors?

All scripts were sourced with a try catch. The first benchmark is
whether or not the code can even be executed. The total number of
executable scripts are summed by model and divided by total evaluations.
*For example, deepseek provided two solutions for problems, this was not
part of the prompt, but both were evaluated, and therefore the Deepseek
execution is based on 10 evals, instead of 5.*

| company   | model                               | execution_ratio |
|-----------|-------------------------------------|-----------------|
| anthropic | claude-opus-4-5-20251101            | 0.6             |
| google    | gemini-3-pro                        | 0.6             |
| deepseek  | deepseek-v3.2-exp-thinking          | 0.5             |
| meta      | llama-4-maverick-03-26-experimental | 0.4             |
| moonshot  | kimi-k2-0905-preview                | 0.4             |
| openai    | gpt-5.1-high                        | 0.4             |
| zai       | glm-4.6                             | 0.4             |
| alibaba   | qwen3-max-preview                   | 0.2             |
| mistral   | mistral-medium-2508                 | 0.2             |
| nous      | hermes-4-405B-reasoning             | 0.2             |
| xai       | grok-4.1                            | 0.2             |
| baidu     | ernie-5.0-preview-1022              | 0.0             |

#### TODO

- 2)  measure correctness for ties (claude-opus-4-5-20251101 vs
      gemini-3-pro)  
- 3)  measure style for ties

## collaboration

- PRs for additional LLMs added to /bench are welcome!
