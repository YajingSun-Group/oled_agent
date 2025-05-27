import os
import json

def doi_encode(doi):
    return doi.replace('/', '%2F')
def doi_decode(doi):
    return doi.replace('%2F', '/')

def read_md(doi, pdf_split_dir="/home/qianzhang/MyProject/deepseek/000-final/pdf_split"):
    encoded_doi = doi_encode(doi)
    md_path = os.path.join(pdf_split_dir, encoded_doi, "output", f"{encoded_doi}.md")
    
    if not os.path.exists(md_path):
        raise FileNotFoundError(f"MD file not found: {md_path}")
        
    with open(md_path, 'r') as f:
        md_content = f.read()
    return md_content

def adjust_if_have_materils(doi):
    encoded_doi = doi_encode(doi)
    extract_json = os.path.join("/home/qianzhang/MyProject/deepseek/000-final/extract_info", encoded_doi, f"{encoded_doi}.json")
    
    # 检查文件是否存在
    if not os.path.exists(extract_json):
        print(f"extract_json not found: {extract_json}")
        return False
        
    # 读取JSON文件
    with open(extract_json, 'r') as f:
        extract_dict = json.load(f)
    
    # 检查是否存在materials键且不为空，并且至少有一个有效材料名称
    has_valid_materials = False
    if 'materials' in extract_dict and extract_dict['materials']:
        for material in extract_dict['materials']:
            if (material.get('emitter_name_full') is not None or 
                material.get('emitter_name_abbreviation') is not None):
                has_valid_materials = True
                break
    
    if not has_valid_materials:
        return False
        
    # 检查devices相关信息
    has_valid_devices = False
    if 'devices' in extract_dict and extract_dict['devices']:
        for device in extract_dict['devices']:
            if device.get('maximum_EQE') is not None:
                has_valid_devices = True
                break
    
    if not has_valid_devices:
        return False
    
    return True
    



