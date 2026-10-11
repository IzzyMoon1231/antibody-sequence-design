"""Synthetic eligibility contracts only; no historical arithmetic or inference."""
import copy
import csv
import importlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import statistics_contract as sc
import run_reproducibility_checks as runner


def sha(text):
    return sc.core.digest_bytes(text.encode())


class ContractFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        code = self.root/'extract.py'; code.write_text('# synthetic extraction identity; never executed\n')
        self.extraction = sc.core.file_sha256(code)
        self.run = {'run_instance_id':'synthetic-run-1','input_bundle_digest':sha('inputs'),
                    'model_id':'synthetic-only','settings':'same-protocol'}
        configuration = sc.metric_configuration(self.run)
        records, eligibility, measurements, endpoints, slots = [], [], [], [], []
        for case in ('9NFU','9NH7'):
            endpoint_id=case+'-H1-plddt'
            endpoint=dict(endpoint_id=endpoint_id,case=case,metric='plddt',region='H1',kind='prediction',
                unit='confidence_0_100',direction='higher',bounds={'lower':0,'upper':100},
                metric_policy_id=sha(case+'policy'),extraction_code_id=self.extraction,
                reference_sha256=sha(case+'reference'),residue_map_sha256=sha('map'),
                definition={'atom_selection':'CA','alignment':'not_applicable','domain':'VH',
                            'mask_source':'current_design_mask','score_definition':'declared synthetic regional confidence'},
                estimands={'primary':'draw_weighted','sensitivity':'unique_sequence_support','unique_support_rule':'retain_all_evaluation_members',
                           'unique_support_weight':'equal_per_sequence_within_case_method',
                           'within_sequence_rule':'mean_over_distinct_prediction_ids'},
                effect={'statistic':'mean_difference','location':'ensemble_minus_baseline','sign':'raw_metric_units'},
                uncertainty={'method':'descriptive_only','unit':'none','hierarchy':'none'},
                replication={'rule':'single_per_draw','assignment':'per_draw','configuration_id':configuration},
                missingness={'rule':'block','estimand':'full_planned_population'},pairing={'mode':'unpaired'})
            endpoints.append(endpoint)
            for method in ('baseline','ensemble'):
                did=sha(case+method+'draw1');pid=sha(did+'prediction1')
                record=dict(design_id=did,record_key=f'{case}/{method}/1',target=case,method=method,
                            source_sequence_sha256=sha(case+method+'sequence'))
                records.append(record)
                observed={**record,'prediction_id':pid,'prediction_replicate_id':'1',
                    'output_pdb_sha256':sha(pid+'bytes'),'reference_sha256':endpoint['reference_sha256'],
                    'residue_map_sha256':sha('map'),'run_contract_sha256':sc.policy_digest(self.run),
                    'prediction_presence_status':'PASS','sequence_validation_status':'PASS',
                    'run_provenance_status':'VERIFIED'}
                eligibility.append(observed)
                measurement={k:observed[k] for k in sc.core.IDENTITY_FIELDS}
                measurement.update(design_id=did,prediction_id=pid,measurement_id=sha(pid+'measurement'),
                    endpoint_id=endpoint_id,metric='plddt',region='H1',unit='confidence_0_100',
                    metric_policy_id=endpoint['metric_policy_id'],extraction_code_id=self.extraction,
                    extraction_run_id='synthetic-extraction-1',extraction_source_path='extract.py',value=80.0,
                    evidence_class='current',source_path='measurements.csv',value_column='value',
                    source_columns={k:k for k in sc.SOURCE_ID_FIELDS})
                measurements.append(measurement)
                slots.append(dict(endpoint_id=endpoint_id,design_id=did,prediction_owner_design_id=did,
                    prediction_id=pid,replicate_id='1',configuration_id=configuration))
        self.contexts=[{'validated':dict(generation_code_commit='a'*40,generation_code_status='VERIFIED',
            mapping_authority_status='VERIFIED',records=records,eligibility=eligibility),'contract':self.run}]
        self.plan=dict(schema_version=sc.PLAN_SCHEMA,analysis_scope='fixed_cases',cases=['9NFU','9NH7'],
            methods=['baseline','ensemble'],design_ids=[r['design_id'] for r in records],endpoints=endpoints,
            evaluation_members=slots,generation_interpretation='implemented_protocol_contrast',analysis_timing='retrospective_exploratory',
            hypotheses=[dict(hypothesis_id=e['endpoint_id']+'-EminusB',case=e['case'],endpoint_id=e['endpoint_id'],
                             comparison=['ensemble','baseline'],family_id='declared-family') for e in endpoints],
            families=[dict(family_id='declared-family',hypothesis_ids=[e['endpoint_id']+'-EminusB' for e in endpoints],
                           correction='descriptive_only',missingness='block')])
        self.measurements=measurements
        self.write_sources()

    def write_sources(self):
        fields=sorted(sc.SOURCE_ID_FIELDS)+['value']
        raw=[]
        for m in self.measurements:
            raw.append({**{k:str(m[k]) for k in sc.SOURCE_ID_FIELDS},'value':'' if m['value'] is None else str(m['value'])})
        path=self.root/'measurements.csv'
        with path.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(raw)
        for i,(m,row) in enumerate(zip(self.measurements,raw),1):
            m.update(source_file_sha256=sc.core.file_sha256(path),source_row_index=i,source_row_sha256=sc.policy_digest(row))

    def approval(self):
        digest=sc.policy_digest(self.plan)
        return digest,dict(schema_version='uirj-plan-approval-v1',plan_sha256=digest,
            decision_id='SYNTHETIC-ONLY',decision_source='synthetic fixture, not a study approval',
            decision_authority='test fixture',decision_timestamp='fixture timestamp')

    def assess(self,approve=True):
        digest,approval=self.approval()
        return sc.assess_validated_context(self.contexts,self.plan,approval if approve else None,
                                           digest if approve else None,self.measurements,self.root)

    def semantics(self, report):
        return {k:v for k,v in report.items() if k not in {'analysis_plan_sha256','approval_binding'}}

    def codes(self,report):return {p['reason_code'] for p in report['reasons']}

    def test_valid_contract_not_scientific_proof(self):
        result=self.assess()
        self.assertEqual(result['state'],'VERIFIED',result['reasons'])
        self.assertEqual(result['design_count'],4);self.assertEqual(len(result['cells']),4)
        self.assertFalse(result['statistics_computed']);self.assertFalse(result['scientific_validity_established'])
        self.assertEqual(result['assurance'],'CONTRACT_ONLY_NOT_SCIENTIFIC_VALIDITY')

    def test_explicit_approval_and_generation_gates(self):
        self.assertIn('SCIENTIFIC_APPROVAL_NOT_SUPPLIED',self.codes(self.assess(False)))
        self.contexts[0]['validated']['generation_code_status']='NOT_ESTABLISHED'
        self.contexts[0]['validated']['generation_code_commit']=None
        report=self.assess()
        self.assertIn('GENERATION_SOURCE_RUN_NOT_ESTABLISHED',self.codes(report))
        self.assertTrue(all(c['scientific_state']!='VERIFIED' for c in report['cells']))
        self.assertTrue(all(c['technical_state']=='VERIFIED' for c in report['cells']))

    def test_nominal_pairing_and_incompatible_bootstrap(self):
        endpoint=self.plan['endpoints'][0]
        endpoint['pairing']={'mode':'paired','basis':'same_seed'}
        self.assertIn('NOMINAL_ID_OR_SEED_PAIRING_REJECTED',self.codes(self.assess()))
        endpoint['pairing']={'mode':'unpaired','pairs':[]}
        self.assertIn('UNPAIRED_PLAN_HAS_PAIRING',self.codes(self.assess()))
        endpoint['pairing']={'mode':'paired','basis':'prespecified_coupling'}
        endpoint['uncertainty']={'method':'bootstrap','unit':'draw','hierarchy':'independent_draws',
                                'resamples':10,'seed':1,'confidence':0.95,'interval_algorithm':'percentile','quantile_rule':'linear'}
        self.assertIn('PAIRED_UNCERTAINTY_HIERARCHY_UNSUPPORTED',self.codes(self.assess()))

    def test_complete_coupling_permits_distinct_sequences(self):
        endpoint=self.plan['endpoints'][0]
        draws={r['method']:r for r in self.contexts[0]['validated']['records'] if r['target']==endpoint['case']}
        pair={'baseline_design_id':draws['baseline']['design_id'],'ensemble_design_id':draws['ensemble']['design_id'],'coupling_id':'real-fixture-coupling'}
        receipt={'coupling_protocol_id':'fixture-protocol','pairs':[pair]}
        path=self.root/'coupling.json';path.write_text(json.dumps(receipt))
        endpoint['pairing']={**receipt,'mode':'paired','basis':'prespecified_coupling','evidence_path':'coupling.json','evidence_sha256':sc.core.file_sha256(path)}
        self.assertNotEqual(draws['baseline']['source_sequence_sha256'],draws['ensemble']['source_sequence_sha256'])
        self.assertEqual(self.assess()['state'],'VERIFIED')

    def test_missing_method_metric_replicate_and_family_preserved(self):
        self.measurements.pop()
        report=self.assess()
        self.assertEqual(report['design_count'],4);self.assertEqual(len(report['cells']),4)
        self.assertIn('MEASUREMENT_MISSING',self.codes(report))
        self.assertEqual(len(report['families'][0]['planned_hypothesis_ids']),2)
        self.assertTrue(report['families'][0]['unavailable_hypothesis_ids'])
        self.plan['evaluation_members'].pop()
        self.assertIn('PLANNED_REPLICATE_MISSING',self.codes(self.assess()))
        self.plan['families'][0]['hypothesis_ids'].pop()
        self.assertIn('HYPOTHESIS_FAMILY_MEMBERSHIP_MISMATCH',self.codes(self.assess()))

    def test_source_row_same_value_wrong_identity_rejected(self):
        a,b=self.measurements[:2]
        for k in ('source_row_index','source_row_sha256'):a[k]=b[k]
        self.assertEqual(a['value'],b['value'])
        self.assertIn('MEASUREMENT_SOURCE_MISMATCH',self.codes(self.assess()))

    def test_source_core_hash_mask_region_and_unit_mutations(self):
        original=copy.deepcopy(self.measurements)
        for field,value in [('source_sequence_sha256',sha('wrong')),('reference_sha256',sha('other reference')),
                            ('residue_map_sha256',sha('wrong map')),('region','other-loop'),('unit','fraction')]:
            with self.subTest(field=field):
                self.measurements=copy.deepcopy(original);self.measurements[0][field]=value;self.write_sources()
                self.assertNotEqual(self.assess()['state'],'VERIFIED')

    def test_historical_null_nonfinite_and_inherent_range(self):
        self.measurements[0]['evidence_class']='historical'
        self.assertIn('HISTORICAL_MEASUREMENT_NOT_PROMOTED',self.codes(self.assess()))
        self.measurements[0]['evidence_class']='current';self.measurements[0]['value']=None
        self.measurements[0]['missing_reason']='extraction unavailable';self.write_sources()
        self.assertIn('MEASUREMENT_NULL',self.codes(self.assess()))
        self.measurements[0]['value']=float('nan');self.write_sources()
        self.assertNotEqual(self.assess()['state'],'VERIFIED')
        self.measurements[0]['value']=101.0;self.write_sources()
        self.plan['endpoints'][0]['bounds']={'lower':None,'upper':None}
        self.assertIn('INHERENT_UNIT_BOUNDS_UNRESOLVED',self.codes(self.assess()))
        self.plan['endpoints'][0]['bounds']={'lower':-10,'upper':1000}
        self.assertIn('PLAN_BOUNDS_EXCEED_UNIT_DOMAIN',self.codes(self.assess()))

    def test_duplicate_measurements_quarantine_all_order_invariant(self):
        bad=copy.deepcopy(self.measurements[0]);bad['value']=90
        self.measurements.append(bad)
        first=self.assess();self.measurements.reverse();second=self.assess()
        self.assertEqual(first,second)
        self.assertIn('DUPLICATE_OR_CONFLICTING_MEASUREMENT',self.codes(first))
        self.assertIn(bad['measurement_id'],first['quarantined_measurement_ids'])

    def test_duplicate_contexts_and_conflicting_run_instances(self):
        duplicate=copy.deepcopy(self.contexts[0]);self.contexts.append(duplicate)
        self.assertIn('DUPLICATE_PREDICTION_CONTEXT',self.codes(self.assess()))
        first=self.assess();self.contexts.reverse();self.assertEqual(self.semantics(first),self.semantics(self.assess()))
        self.contexts[1]['contract']['settings']='different-model'
        self.assertIn('CONFLICTING_RUN_INSTANCE_DEFINITION',self.codes(self.assess()))

    def test_permutations_and_unplanned_rows_propagate(self):
        first=self.assess()
        self.measurements.reverse();self.plan['evaluation_members'].reverse()
        self.plan['endpoints'].reverse();self.plan['hypotheses'].reverse()
        self.contexts[0]['validated']['records'].reverse();self.contexts[0]['validated']['eligibility'].reverse()
        self.assertEqual(self.semantics(first),self.semantics(self.assess()))
        extra=copy.deepcopy(self.measurements[0]);extra['endpoint_id']='not-in-plan';extra['measurement_id']=sha('extra')
        self.measurements.append(extra);self.write_sources()
        report=self.assess()
        self.assertIn('UNPLANNED_MEASUREMENT_RETAINED',self.codes(report))
        self.assertTrue(all(c['scientific_state']!='VERIFIED' for c in report['cells']))
        self.assertTrue(report['families'][0]['unavailable_hypothesis_ids'])

    def test_conflicting_endpoint_permutation_and_inference_gate(self):
        duplicate=copy.deepcopy(self.plan['endpoints'][0]);duplicate['unit']='nonsense'
        self.plan['endpoints'].append(duplicate)
        first=self.assess();self.plan['endpoints'].reverse();self.assertEqual(self.semantics(first),self.semantics(self.assess()))
        self.assertIn('DUPLICATE_ENDPOINT_ID',self.codes(first))
        self.plan['endpoints']=[e for e in self.plan['endpoints'] if e['unit']!='nonsense']
        self.plan['families'][0]['correction']='holm'
        self.assertIn('INFERENTIAL_HYPOTHESIS_CONTRACT_UNRESOLVED',self.codes(self.assess()))

    def test_explicit_effect_unique_aggregation_and_nested_shapes(self):
        endpoint=self.plan['endpoints'][0]
        original=copy.deepcopy(endpoint)
        endpoint['effect'].pop('statistic')
        self.assertIn('EFFECT_DEFINITION_UNSUPPORTED_OR_UNRESOLVED',self.codes(self.assess()))
        self.plan['endpoints'][0]=copy.deepcopy(original)
        self.plan['endpoints'][0]['estimands'].pop('within_sequence_rule')
        self.assertIn('ESTIMAND_POLICY_UNSUPPORTED_OR_UNRESOLVED',self.codes(self.assess()))
        for field in ('definition','estimands','uncertainty','replication','pairing'):
            with self.subTest(field=field):
                self.plan['endpoints'][0]=copy.deepcopy(original)
                self.plan['endpoints'][0][field]=None
                self.assertIn('MALFORMED_ANALYSIS_PLAN',self.codes(self.assess()))

    def test_weak_mean_permutation_remains_unsupported(self):
        self.plan['families'][0].update(correction='holm',alpha=0.05)
        for hypothesis in self.plan['hypotheses']:
            hypothesis.update(test='permutation_mean_difference',null='zero_mean_difference',
                              sidedness='two_sided',sampling_assumptions='arbitrary unequal distributions')
        self.assertIn('INFERENTIAL_HYPOTHESIS_CONTRACT_UNRESOLVED',self.codes(self.assess()))

    def test_repeated_draw_shared_prediction_keeps_draw_weights(self):
        original=self.contexts[0]['validated']['records'][0]
        duplicate=dict(original,design_id=sha('separate-generated-draw'),record_key='9NFU/baseline/2')
        self.contexts[0]['validated']['records'].append(duplicate)
        self.plan['design_ids'].append(duplicate['design_id'])
        member=copy.deepcopy(self.plan['evaluation_members'][0]);member['design_id']=duplicate['design_id']
        self.plan['evaluation_members'].append(member)
        self.plan['endpoints'][0]['replication'].update(rule='shared_deterministic',assignment='shared_sequence')
        report=self.assess()
        self.assertEqual(report['state'],'VERIFIED',report['reasons'])
        self.assertEqual(report['design_count'],5);self.assertEqual(report['unique_sequence_count'],4)
        self.assertEqual(report['prediction_identity_count'],4)
        self.assertEqual(len(report['estimand_memberships']['draw_weighted']),5)
        support=[s for s in report['estimand_memberships']['unique_sequence_support'] if s['source_sequence_sha256']==original['source_sequence_sha256']][0]
        self.assertEqual(len(support['design_ids']),2);self.assertEqual(len(support['prediction_ids']),1)
        self.plan['endpoints'][0]['replication'].update(rule='single_per_draw',assignment='per_draw')
        self.assertIn('REUSE_MISLABELED_AS_INDEPENDENT_FOLD',self.codes(self.assess()))

    def test_independent_prediction_replicates_kept_and_configuration_checked(self):
        second_run=dict(self.run,run_instance_id='synthetic-run-2',input_bundle_digest=sha('replicate2 inputs'))
        observed=copy.deepcopy(self.contexts[0]['validated']['eligibility'][0])
        observed.update(prediction_id=sha('independent-fold-2'),prediction_replicate_id='2',
                        run_contract_sha256=sc.policy_digest(second_run),output_pdb_sha256=sha('fold2 bytes'))
        context={'contract':second_run,'validated':{**self.contexts[0]['validated'],
            'records':copy.deepcopy(self.contexts[0]['validated']['records']),'eligibility':[observed]}}
        self.contexts.append(context)
        measurement=copy.deepcopy(self.measurements[0]);measurement.update({k:observed[k] for k in sc.core.IDENTITY_FIELDS})
        measurement.update(prediction_id=observed['prediction_id'],measurement_id=sha('independent-measurement2'),value=75.0)
        self.measurements.append(measurement);self.write_sources()
        member=copy.deepcopy(self.plan['evaluation_members'][0]);member.update(prediction_id=observed['prediction_id'],replicate_id='2')
        self.plan['evaluation_members'].append(member)
        self.plan['endpoints'][0]['replication']['rule']='mean_within_draw'
        report=self.assess()
        self.assertEqual(report['state'],'VERIFIED',report['reasons'])
        self.assertEqual(report['design_count'],4);self.assertEqual(report['prediction_identity_count'],5)
        cell=[c for c in report['cells'] if c['design_id']==observed['design_id']][0]
        self.assertEqual(len(cell['retained_members']),2)
        self.assertEqual(len({m['prediction_id'] for m in cell['retained_members']}),2)
        self.contexts[1]['contract']['settings']='different-model-config'
        self.assertIn('REPLICATE_CONFIGURATION_MISMATCH',self.codes(self.assess()))

    def test_sequence_score_remains_unsupported(self):
        self.plan['endpoints'][0]['kind']='sequence_score'
        self.assertIn('SEQUENCE_SCORE_PROVENANCE_UNSUPPORTED',self.codes(self.assess()))

    def test_shared_validator_rejection_and_no_plan_coverage(self):
        with mock.patch.object(sc.core,'validate_prediction_bundle',side_effect=ValueError('actual shared rejection')):
            with self.assertRaisesRegex(ValueError,'actual shared rejection'):
                sc.evaluate_eligibility(['fixture'],self.root,self.root,self.root)
        with mock.patch.object(sc.core,'validate_prediction_bundle',return_value=self.contexts[0]['validated']) as shared:
            result=sc.evaluate_eligibility(['fixture'],self.root,self.root,self.root)
        shared.assert_called_once_with('fixture',self.root,self.root,self.root,None)
        self.assertIn('ANALYSIS_PLAN_MISSING',self.codes(result))
        self.assertEqual(result['prediction_coverage'][0]['prediction_presence_counts'],{'PASS':4})

    def test_import_modes_cannot_accidentally_compute(self):
        with mock.patch.object(runner.subprocess,'run',side_effect=AssertionError('legacy computation forbidden')):
            importlib.reload(runner)
            with self.assertRaises(SystemExit):runner.main([])
            report=self.assess()
            with mock.patch.object(sc,'evaluate_eligibility',return_value=report),mock.patch('sys.stdout',new_callable=io.StringIO):
                self.assertEqual(runner.main(['--mode','eligibility','--manifest','fixture',
                    '--input-root','input','--upstream-root','upstream','--output-root','output']),0)
        fake=mock.Mock(returncode=0,stdout='synthetic check',stderr='')
        with mock.patch.object(runner.subprocess,'run',return_value=fake) as run,mock.patch('sys.stdout',new_callable=io.StringIO):
            self.assertEqual(runner.main(['--mode','legacy-arithmetic']),0)
        self.assertEqual(run.call_count,4)
        with mock.patch.object(runner.subprocess,'run',side_effect=AssertionError('forbidden')):
            with self.assertRaises(SystemExit):runner.main(['--mode','legacy-arithmetic','--manifest','fixture'])


if __name__=='__main__':
    unittest.main()
