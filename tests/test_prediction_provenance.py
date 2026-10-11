"""Bounded integrity fixtures. No inference, scientific metrics or downloads."""
import copy
import csv
import hashlib
import importlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import validate_esmfold_predictions as v
import build_prediction_manifest as b
import run_esmfold_predictions as runner


def ca_pdb(sequence='G', model=None, insertion=' ', duplicate=False, absent=False):
    inverse = {a: name for name, a in v.AA3.items()}
    lines = [f'MODEL     {model:4d}\n'] if model is not None else []
    for i, aa in enumerate(sequence, 1):
        atom = 'N' if absent else 'CA'
        line = f'ATOM  {i:5d} {atom:^4s} {inverse[aa]:3s} A{i:4d}{insertion}   {1.0:8.3f}{2.0:8.3f}{3.0:8.3f}{1.0:6.2f}{80.0:6.2f}           C  \n'
        lines.append(line)
        if duplicate:
            lines.append(line)
    lines.append('ENDMDL\n' if model is not None else 'END\n')
    return ''.join(lines)


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.local = self.root/'local'; self.upstream = self.root/'upstream'
        self.out = self.root/'outputs'
        self.local.mkdir(); self.upstream.mkdir()
        self.sequence = 'G'*115
        content = f'>T=0.1, sample=1, score=0\n{self.sequence}\n'
        for root in (self.local, self.upstream):
            (root/'input.fa').write_text(content)
        sha = v.file_sha256(self.local/'input.fa')
        row = dict(schema_version=v.SCHEMA, target='9NH7', method='baseline', index='1',
            sequence=self.sequence, pdb_path='draft/9NH7/baseline/1/replicate-1.pdb',
            record_key='9NH7/baseline/1', source_header='T=0.1, sample=1, score=0',
            source_design_id='1', source_record_ordinal='1', source_sequence_sha256=v.digest_bytes(self.sequence.encode()),
            expected_length='115', source_repo=b.REPO2, source_commit='a'*40, source_path='input.fa',
            source_file_sha256=sha, upstream_repo=b.REPO1, upstream_commit='b'*40,
            upstream_path='input.fa', upstream_file_sha256=sha, reference_id='9NH7:'+'c'*64,
            reference_sha256='c'*64, residue_map_sha256='d'*64, mapping_authority_status=v.AUTHORITY,
            planned_run_id='', prediction_replicate_id='1', prediction_id='')
        row['design_id'] = v.design_identity(row)
        self.rows = [row]
        self.contract = dict(schema_version='uirj-run-v1', model_id='facebook/esmfold_v1',
            model_revision='e'*40, assets={'config.json':'f'*64}, serializer_sha256=v.file_sha256(runner.__file__),
            runtime=runner.observed_runtime_versions(), device='cuda',
            precision={'esm':'float16','trunk':'float32'}, inference_settings={'num_recycles':4,'chunk_size':64},
            seed=0, deterministic_algorithms=True, input_bundle_digest=v.input_bundle_digest(self.rows),
            run_instance_id='fixture-run')

    def plan(self):
        row=self.rows[0]
        run=v.validate_run_contract(self.contract,self.rows)
        row['planned_run_id']=run; row['prediction_id']=v.prediction_identity(row)
        row['pdb_path']=row['pdb_path'].replace('draft/',run+'/',1)
        return row

    def validate(self):
        return v.validate_input_records(self.rows,self.local,self.upstream,self.out)

    def test_source_and_duplicate_sequence_keep_design_identity(self):
        self.validate()
        row=copy.deepcopy(self.rows[0]);row['index']='2';row['source_record_ordinal']='2'
        row['source_design_id']='2';row['record_key']='9NH7/baseline/2'
        row['source_header']='T=0.1, sample=2, score=0'
        row['pdb_path']='draft/9NH7/baseline/2/replicate-1.pdb'
        for root in (self.local,self.upstream):
            with (root/'input.fa').open('a') as h:h.write(f">{row['source_header']}\n{self.sequence}\n")
        self.rows.append(row)
        for r in self.rows:
            r['source_file_sha256']=r['upstream_file_sha256']=v.file_sha256(self.local/'input.fa')
            r['design_id']=v.design_identity(r)
        self.validate()
        self.assertNotEqual(self.rows[0]['design_id'],row['design_id'])
        self.assertEqual(self.rows[0]['source_sequence_sha256'],row['source_sequence_sha256'])
        self.assertEqual(v.input_bundle_digest(self.rows),v.input_bundle_digest(list(reversed(self.rows))))

    def test_strict_wrapper_rejects_subset_before_map(self):
        manifest=self.root/'manifest.csv'
        with manifest.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=b.MANIFEST_FIELDS)
            writer.writeheader();writer.writerows(self.rows)
        with self.assertRaisesRegex(ValueError,'Incomplete frozen'):
            v.load_input_manifest(manifest,self.local,self.upstream,self.out)

    def test_input_mutations_rejected(self):
        mutations=[('record_key','wrong'),('source_header','sample=1 sample=1'),
          ('source_sequence_sha256','0'*64),('sequence','X'*115),('expected_length','114'),
          ('source_record_ordinal','2'),('source_design_id','0'),('pdb_path','../escape'),
          ('source_path','/absolute'),('mapping_authority_status','APPROVED'),
          ('prediction_id','invented')]
        original=copy.deepcopy(self.rows[0])
        for field,value in mutations:
            with self.subTest(field=field):
                self.rows[0]=dict(original,**{field:value})
                with self.assertRaises(ValueError):self.validate()
        self.rows=[original,copy.deepcopy(original)]
        with self.assertRaises(ValueError):self.validate()

    def test_changed_source_and_symlink_escape_rejected(self):
        (self.local/'input.fa').write_text('>sample=1\nA\n')
        with self.assertRaises(ValueError):self.validate()
        (self.local/'escape').symlink_to(self.root)
        with self.assertRaises(ValueError):v.rooted_path(self.local,'escape/out.pdb')
        for path in ('../out','/out','C:/out','a/../b','a\\b'):
            with self.assertRaises(ValueError):v.rooted_path(self.out,path)

    def test_fasta_invalid_symbols_headers_and_whitespace(self):
        p=self.root/'bad.fa'
        for content in ('>\nG\n','G\n','>one\nG G\n','>one\ng\n','>one\nX\n'):
            p.write_text(content)
            with self.assertRaises(ValueError):v.read_fasta_records(p)

    def test_pdb_model_insertion_absent_ca_and_duplicate(self):
        p=self.root/'reference.pdb';p.write_text(ca_pdb(model=7,insertion='A'))
        record=v.parse_pdb_residues(p)[0]
        self.assertEqual((record['model'],record['insertion_code']),(7,'A'))
        p.write_text(ca_pdb(absent=True));self.assertFalse(v.parse_pdb_residues(p)[0]['ca_present'])
        p.write_text(ca_pdb(duplicate=True))
        with self.assertRaises(ValueError):v.parse_pdb_residues(p)

    def test_map_keeps_gaps_and_atom_seqres_discrepancy(self):
        p=self.root/'ref.pdb'
        # Small coordinate fixture: SEQRES A/G, observed atom G at position1, position2 absent.
        text=ca_pdb('G').replace(' A   1',' C   1')
        p.write_text('SEQRES   1 C    2  ALA GLY\n'+text)
        rows=b.map_reference(p,'9NFU','C',2,1,0,{0},{1},[('H',0,None)],[('L',0,None)],'a'*40)
        self.assertEqual(rows[0]['reference_atom_seqres_agree'],'false')
        self.assertEqual(rows[1]['fill_required'],'true')
        self.assertEqual(rows[1]['reference_atom_aa'],'')
        self.assertEqual((rows[0]['current_domain'],rows[0]['historical_domain']),('H','L'))
        with self.assertRaises(ValueError):b.domain_at([('H',0,2),('L',0,2)],0,2)
        p.write_text('SEQRES   1 C    2  ALA GLY\n'+ca_pdb('G',insertion='A').replace(' A   1',' C   1'))
        with self.assertRaises(ValueError):b.map_reference(p,'9NFU','C',2,1,0,{0},{1},[('H',0,None)],[('L',0,None)],'a'*40)

    def test_prediction_replicates_and_run_registry(self):
        row=self.plan();one=row['prediction_id']
        row['prediction_replicate_id']='2'
        self.assertNotEqual(one,v.prediction_identity(row))
        with self.assertRaises(ValueError):v.validate_run_contract(self.contract,self.rows)
        row['prediction_replicate_id']='1'
        v.validate_run_registry([self.contract],self.rows)
        with self.assertRaises(ValueError):v.validate_run_registry([self.contract,self.contract],self.rows)

    def test_receipt_storage_and_stale_run(self):
        row=self.plan();path=v.rooted_path(self.out,row['pdb_path'])
        context={'row':row,'contract':self.contract,'rows':self.rows}
        runner.save_prediction(ca_pdb(self.sequence),path,self.sequence,context)
        receipt=json.loads(path.with_name(path.name+'.provenance.json').read_text())
        v.validate_receipt(receipt,row,self.contract,self.rows,path)
        result=v._prediction_eligibility(row,self.out,self.rows,self.contract)
        self.assertEqual(result['run_provenance_status'],'VERIFIED')
        self.assertEqual(result['final_metric_eligibility'],'NOT_ASSESSED')
        for field in ('record_key','prediction_id','source_sequence_sha256','output_pdb_sha256','run_contract_sha256'):
            bad=dict(receipt,**{field:'wrong'})
            with self.assertRaises(ValueError):v.validate_receipt(bad,row,self.contract,self.rows,path)
        with self.assertRaises(FileExistsError):runner.save_prediction(ca_pdb(self.sequence),path,self.sequence,context)
        path.with_name(path.name+'.provenance.json').unlink()
        self.assertEqual(v._prediction_eligibility(row,self.out,self.rows,self.contract)['run_provenance_status'],'RUN_UNVERIFIED')
        with self.assertRaises(ValueError):v.prediction_eligibility(row,self.out,self.rows,self.contract)

    def test_fake_runtime_receipt_rejected_before_write(self):
        for field, value in [('runtime', {'python':'test','torch':'test','transformers':'test'}),
                             ('serializer_sha256', '1'*64)]:
            with self.subTest(field=field):
                original=copy.deepcopy(self.contract)
                self.contract[field]=value
                row=self.plan();path=v.rooted_path(self.out,row['pdb_path'])
                with self.assertRaisesRegex(ValueError, 'actual runtime/serializer'):
                    runner.save_prediction(ca_pdb(self.sequence),path,self.sequence,
                        {'row':row,'contract':self.contract,'rows':self.rows})
                self.assertFalse(path.exists())
                self.contract=original

    def test_source_tuple_json_list_semantic_equivalence(self):
        literal={'C':[(None,0,3),('H',3,124),(None,124,143),('L',143,None)]}
        serialized=json.loads(json.dumps(literal))
        self.assertNotEqual(literal,serialized)
        self.assertEqual(v.canonical_bytes(literal),v.canonical_bytes(serialized))
        serialized['C'][1][1]=9
        self.assertNotEqual(v.canonical_bytes(literal),v.canonical_bytes(serialized))

    def test_duplicate_ca_never_strict_pass(self):
        row=self.plan();path=v.rooted_path(self.out,row['pdb_path']);path.parent.mkdir(parents=True)
        path.write_text(ca_pdb(self.sequence,duplicate=True))
        path.with_name(path.name+'.complete').write_text('COMPLETE\n')
        self.assertNotEqual(v._prediction_eligibility(row,self.out,self.rows,self.contract)['sequence_validation_status'],'PASS')

    def test_interrupted_receipt_keeps_incomplete_file(self):
        row=self.plan();path=v.rooted_path(self.out,row['pdb_path'])
        with mock.patch.object(v,'validate_receipt',side_effect=ValueError('injected interruption')):
            with self.assertRaises(ValueError):runner.save_prediction(ca_pdb(self.sequence),path,self.sequence,{'row':row,'contract':self.contract,'rows':self.rows})
        self.assertTrue(path.is_file());self.assertFalse(path.with_name(path.name+'.complete').exists())
        self.assertEqual(runner.check_prediction(path,self.sequence),'INCOMPLETE_NO_MARKER')

    def test_measurement_ambiguity_and_unverified_identity(self):
        row=self.plan();path=v.rooted_path(self.out,row['pdb_path'])
        runner.save_prediction(ca_pdb(self.sequence),path,self.sequence,{'row':row,'contract':self.contract,'rows':self.rows})
        identity=v._prediction_eligibility(row,self.out,self.rows,self.contract)
        measurement={k:identity[k] for k in v.IDENTITY_FIELDS}
        measurement.update(design_id=row['design_id'],prediction_id=row['prediction_id'],measurement_id='m1',metric='RMSD',region='H1',metric_policy_id='unapproved-test',extraction_code_id='test-only')
        v.validate_measurement_joins([measurement],[identity])
        with self.assertRaises(ValueError):v.validate_measurement_joins([measurement,measurement],[identity])
        with self.assertRaises(ValueError):v.validate_measurement_joins([measurement],[identity,identity])
        with self.assertRaises(ValueError):v.validate_measurement_joins([dict(measurement,run_contract_sha256='f'*64)],[identity])
        with self.assertRaises(ValueError):v.validate_measurement_joins([measurement],[dict(identity,run_provenance_status='RUN_UNVERIFIED')])

    def test_dry_run_never_model_or_output(self):
        row=self.rows[0]
        with mock.patch.object(runner,'load_manifest',return_value=self.rows), mock.patch.object(runner,'load_model',side_effect=AssertionError('model forbidden')), mock.patch.object(sys,'argv',['runner','--manifest','fixture','--output-root',str(self.out),'--dry-run']):
            runner.main()
        self.assertFalse(self.out.exists())
        with mock.patch.object(runner,'load_manifest',return_value=self.rows), mock.patch.object(runner,'load_model',side_effect=AssertionError('model forbidden')), mock.patch.object(sys,'argv',['runner','--manifest','fixture','--output-root',str(self.out)]):
            with self.assertRaises(ValueError):runner.main()
        self.assertFalse(self.out.exists())

    def test_local_assets_complete_inventory_and_runtime(self):
        row=self.plan();assets=self.root/'assets';assets.mkdir();(assets/'config.json').write_text('{}')
        self.contract['assets']={'config.json':v.file_sha256(assets/'config.json')}
        self.contract['runtime']={'python':'test','torch':'test','transformers':'test'}
        self.contract['serializer_sha256']=v.file_sha256(runner.__file__)
        self.contract['input_bundle_digest']=v.input_bundle_digest(self.rows)
        row['planned_run_id']=v.validate_run_contract(self.contract,self.rows);row['prediction_id']=v.prediction_identity(row)
        with mock.patch('platform.python_version',return_value='test'),mock.patch('importlib.metadata.version',return_value='test'):
            runner.verify_execution_contract(self.contract,self.rows,assets)
            (assets/'unlisted.json').write_text('{}')
            with self.assertRaises(ValueError):runner.verify_execution_contract(self.contract,self.rows,assets)




class ModelContractTests(unittest.TestCase):
    def test_parameter_observation(self):
        class Parameter:
            def __init__(self, dtype, device='cuda'):
                self.dtype=dtype; self.device=mock.Mock(type=device)
            def is_floating_point(self): return True
        model=mock.Mock()
        model.named_parameters.return_value=[('esm.weight',Parameter(runner.torch.float16)),('trunk.weight',Parameter(runner.torch.float32))]
        runner.verify_model_parameters(model)
        model.named_parameters.return_value=[('trunk.weight',Parameter(runner.torch.float16))]
        with self.assertRaises(ValueError):runner.verify_model_parameters(model)
        model.named_parameters.return_value=[('trunk.weight',Parameter(runner.torch.float32,'cpu'))]
        with self.assertRaises(ValueError):runner.verify_model_parameters(model)

    def test_load_model_uses_local_pinned_assets_and_chunk(self):
        import types
        tokenizer_loader=mock.Mock()
        model_loader=mock.Mock()
        model=mock.Mock()
        model.eval.return_value=model;model.cuda.return_value=model
        model.esm.half.return_value=model.esm
        model_loader.from_pretrained.return_value=model
        fake=types.SimpleNamespace(AutoTokenizer=tokenizer_loader,EsmForProteinFolding=model_loader)
        contract={'model_revision':'e'*40,'inference_settings':{'chunk_size':64,'num_recycles':4}}
        with mock.patch.dict(sys.modules,{'transformers':fake}), mock.patch.object(runner.torch.cuda,'is_available',return_value=True), mock.patch.object(runner,'verify_model_parameters') as verify:
            loaded,tokenizer=runner.load_model(contract,'/explicit/local/assets')
        self.assertIs(loaded,model)
        tokenizer_loader.from_pretrained.assert_called_once_with('/explicit/local/assets',local_files_only=True,revision='e'*40)
        model_loader.from_pretrained.assert_called_once_with('/explicit/local/assets',local_files_only=True,revision='e'*40,low_cpu_mem_usage=True)
        model.set_chunk_size.assert_called_once_with(64)
        verify.assert_called_once_with(model)
        model.esm.half.assert_called_once()

    def test_predict_passes_resolved_recycles_to_model(self):
        token=mock.Mock();token.cuda.return_value=token
        tokenizer=mock.Mock(return_value={'input_ids':token})
        model=mock.Mock(return_value={'fixture':'output'})
        with mock.patch.object(runner,'serialize_esmfold_output',return_value='fixture PDB') as serialize:
            result=runner.predict('G',model,tokenizer,{'chunk_size':64,'num_recycles':4})
        self.assertEqual(result,'fixture PDB')
        tokenizer.assert_called_once_with(['G'],return_tensors='pt',add_special_tokens=False)
        model.assert_called_once_with(input_ids=token,num_recycles=4)
        serialize.assert_called_once_with({'fixture':'output'})

    def test_lock_inspector_legacy(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); dest=root/'p.pdb'; lock=root/'p.pdb.lock'
            lock.write_text('')
            manifest=root/'manifest.csv'
            with manifest.open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=['target','method','index','sequence','pdb_path'])
                writer.writeheader();writer.writerow(dict(target='9NH7',method='baseline',index=0,sequence='G',pdb_path=str(dest)))
            result=runner.inspect_prediction_locks(manifest)
            self.assertEqual(result[0]['prediction_status'],'MISSING')
            self.assertTrue(lock.exists())
            with manifest.open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=['schema_version','target','method','index','sequence','pdb_path'])
                writer.writeheader();writer.writerow(dict(schema_version=v.SCHEMA,target='9NH7',method='baseline',index=1,sequence='G',pdb_path='p.pdb'))
            with self.assertRaisesRegex(ValueError,'Explicit root required'):
                runner.inspect_prediction_locks(manifest)
            strict=runner.inspect_prediction_locks(manifest,output_root=root)
            self.assertEqual(strict[0]['prediction_status'],'MISSING')
            self.assertEqual(strict[0]['lock_path'],str(lock.resolve()))

if __name__ == '__main__':
    unittest.main()
