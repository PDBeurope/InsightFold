from pathlib import Path
import tempfile
import unittest

from spec_evaluator import (
    Candidate,
    EvaluationError,
    ReferenceSpec,
    build_rankings,
    discover_candidates,
    evaluate_directory,
    load_reference,
    normalize_judge_response,
    parse_candidate_filename,
)


def good_response(reference_alignment=4):
    return {
        "scores": {
            "correctness": 4,
            "completeness": 3,
            "specificity": 5,
            "coherence": 4,
            "actionability": 3,
            "reference_alignment": reference_alignment,
        },
        "justification": "Solid but incomplete.",
        "strengths": ["Clear structure"],
        "weaknesses": ["Some missing edge cases"],
        "suggested_improvements": ["Add acceptance criteria"],
        "missing_from_reference": ["Reference-only workflow detail"],
        "extra_or_conflicting_claims": ["No conflicts"],
    }


class FilenameParsingTests(unittest.TestCase):
    def test_parse_explicit_separator_with_hyphens_and_underscores(self):
        project, model = parse_candidate_filename(
            Path("protein-interaction_energy--gemma-31b_v2.md")
        )

        self.assertEqual(project, "protein-interaction_energy")
        self.assertEqual(model, "gemma-31b_v2")

    def test_rejects_missing_separator(self):
        with self.assertRaises(ValueError):
            parse_candidate_filename(Path("spec_nmr_restraints_gemma_31b.md"))


class DiscoveryAndReferenceTests(unittest.TestCase):
    def test_discover_groups_candidates_by_project(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "alpha--model_a.md").write_text("A", encoding="utf-8")
            (root / "alpha--model_b.md").write_text("B", encoding="utf-8")
            (root / "beta--model_a.md").write_text("C", encoding="utf-8")

            grouped = discover_candidates(root)

            self.assertEqual(set(grouped), {"alpha", "beta"})
            self.assertEqual([candidate.model for candidate in grouped["alpha"]], ["model_a", "model_b"])

    def test_discover_rejects_invalid_markdown_names_in_strict_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "alpha--model_a.md").write_text("A", encoding="utf-8")
            (root / "alpha_model_b.md").write_text("B", encoding="utf-8")

            with self.assertRaises(ValueError):
                discover_candidates(root, strict=True)

    def test_reference_lookup_present_and_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "alpha.md").write_text("Reference", encoding="utf-8")

            reference = load_reference(root, "alpha")
            missing = load_reference(root, "beta")

            self.assertIsNotNone(reference)
            self.assertEqual(reference.markdown, "Reference")
            self.assertIsNone(missing)


class ValidationTests(unittest.TestCase):
    def test_normalizes_response_and_computes_weighted_score_with_reference(self):
        scores, overall, normalized = normalize_judge_response(
            good_response(reference_alignment=5),
            used_reference=True,
        )

        self.assertEqual(scores["reference_alignment"], 5)
        self.assertEqual(overall, 4.08)
        self.assertEqual(normalized["strengths"], ["Clear structure"])

    def test_reference_alignment_is_null_without_reference(self):
        response = good_response(reference_alignment=None)
        scores, overall, _ = normalize_judge_response(
            response,
            used_reference=False,
        )

        self.assertIsNone(scores["reference_alignment"])
        self.assertEqual(overall, 3.8)

    def test_rejects_bad_score(self):
        response = good_response()
        response["scores"]["correctness"] = 6

        with self.assertRaises(EvaluationError):
            normalize_judge_response(response, used_reference=True)


class BatchEvaluationTests(unittest.TestCase):
    def test_evaluate_directory_with_mock_judge_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates"
            references = root / "references"
            output = root / "results"
            candidates.mkdir()
            references.mkdir()

            (candidates / "alpha--model_a.md").write_text("Candidate A", encoding="utf-8")
            (candidates / "alpha--model_b.md").write_text("Candidate B", encoding="utf-8")
            (references / "alpha.md").write_text("Reference", encoding="utf-8")

            def judge(candidate: Candidate, reference: ReferenceSpec | None):
                response = good_response(reference_alignment=4)
                if candidate.model == "model_b":
                    response["scores"]["correctness"] = 5
                return response

            result = evaluate_directory(
                candidates,
                output,
                reference_dir=references,
                model="mock-model",
                judge=judge,
            )

            self.assertTrue((output / "evaluations.json").exists())
            self.assertTrue((output / "summary.md").exists())
            self.assertTrue((output / "plots" / "alpha_category_bars.svg").exists())
            self.assertTrue((output / "plots" / "alpha_score_heatmap.svg").exists())
            self.assertEqual(len(result["evaluations"]), 2)
            self.assertIn("plots", result)
            self.assertEqual(result["rankings"]["alpha"][0]["model"], "model_b")

    def test_build_rankings_puts_errors_last(self):
        from spec_evaluator import EvaluationRecord

        records = [
            EvaluationRecord(
                project="alpha",
                model="bad",
                candidate_path="bad.md",
                reference_path=None,
                used_reference=False,
                status="error",
                scores={},
                overall_score=None,
                justification="",
                strengths=[],
                weaknesses=[],
                suggested_improvements=[],
                missing_from_reference=[],
                extra_or_conflicting_claims=[],
                error="failed",
            ),
            EvaluationRecord(
                project="alpha",
                model="good",
                candidate_path="good.md",
                reference_path=None,
                used_reference=False,
                status="success",
                scores={},
                overall_score=4.2,
                justification="",
                strengths=[],
                weaknesses=[],
                suggested_improvements=[],
                missing_from_reference=[],
                extra_or_conflicting_claims=[],
            ),
        ]

        ranking = build_rankings(records)

        self.assertEqual(ranking["alpha"][0]["model"], "good")
        self.assertEqual(ranking["alpha"][1]["model"], "bad")


if __name__ == "__main__":
    unittest.main()
