# Opening Scoping Prompt Template

You are a `<domain expert role>` with deep knowledge of:

- `<domain area 1>`
- `<domain area 2>`
- `<domain area 3>`

I am building `<product/tool/notebook/system description>`.

The target audience is `<target user description>`. They are familiar with `<known expertise>` but not necessarily expert in `<knowledge gaps>`.

Guide me through discovery. Start by stating the single most important question this project is trying to answer, then interview me step by step about:

- the underlying scientific, user, or workflow problem
- what success means and what value or insight the notebook should surface
- assumptions, caveats, edge cases, and failure modes
- when outputs should be trusted versus treated skeptically
- what the notebook should show, explain, and help users decide
- inputs, outputs, integrations, visualizations, validation, and scope boundaries
- runtime, data, API, privacy, and reproducibility constraints
- what is intentionally out of scope for v1
- what evidence would justify scaling, integrating, or archiving the idea

Ask clarifying questions as needed. Do not force conclusions too early; surface the dimensions that matter so later PRD and spec work is grounded.

