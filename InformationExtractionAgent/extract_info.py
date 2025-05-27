import pandas as pd
import os
from openai import OpenAI
import base64
import json
from multiprocessing import Pool
from tqdm import tqdm
from llm import LLMCaller, VisualLLMCaller, get_text_prompt
from utils import read_md, doi_encode, doi_decode
import time

# API密钥列表
API_KEYS = [

]

# 基础配置
BASE_CONFIG = {
    "model": "deepseek-chat",
    "base_url": "https://api.deepseek.com/v1",
}

def get_config_with_key(api_key):
    """获取包含特定API密钥的配置"""
    return {
        **BASE_CONFIG,
        "api_key": api_key
    }

def process_single_doi(args):
    """处理单个DOI的函数"""
    doi, api_key = args
    try:
        if os.path.exists(f"/home/qianzhang/MyProject/deepseek/000-final/extract_info/{doi_encode(doi)}/{doi_encode(doi)}.json"):
            print(f"文件已存在: {doi}.json")
            return doi, True, None

        time.sleep(0.1)  # 添加小延迟以避免请求过快
        
        # 使用指定的API密钥创建LLM实例
        config = get_config_with_key(api_key)
        llm = LLMCaller(
            model=config["model"], 
            api_key=config["api_key"], 
            base_url=config["base_url"],
            # temperature=0.1,  # 添加temperature参数
            # top_p=0.95       # 添加top_p参数
        )
        
        # 读取文本并生成提示词
        doi_text = read_md(doi)
        prompt = get_text_prompt(doi_text)
        
        # 调用LLM
        result = llm.call_llm(prompt, response_json=True, stream=False)
        
        # 解析结果
        result_dict = json.loads(result)
        
        # 保存结果
        save_dir = f"/home/qianzhang/MyProject/deepseek/000-final/extract_info/{doi_encode(doi)}"
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, f"{doi_encode(doi)}.json")
        
        with open(save_path, "w", encoding='utf-8') as f:
            json.dump(result_dict, f, indent=2, ensure_ascii=False)
            
        return doi, True, None
    except Exception as e:
        return doi, False, str(e)

def main():
    # 读取DOI列表
    df = pd.read_csv('/home/qianzhang/MyProject/deepseek/000-final/scripts/files/split_pdfs_info-20250417.csv')
    df = df.sort_values(by='年份', ascending=False)
    
    doi_list = df['DOI'].tolist()
    doi_list = [str(x) for x in doi_list if x is not None][10000:20000]
    
    # 创建参数列表，为每个DOI分配一个API密钥
    args_list = []
    for i, doi in enumerate(doi_list):
        # 循环使用API密钥
        api_key = API_KEYS[i % len(API_KEYS)]
        args_list.append((doi, api_key))
    
    # 创建进程池，进程数等于API密钥数量
    n_processes = len(API_KEYS)
    print(f"使用 {n_processes} 个进程进行并行处理")
    
    # 记录成功和失败的DOI
    success_dois = []
    failed_dois = []
    
    # 使用进程池处理
    with Pool(n_processes) as pool:
        results = list(tqdm(
            pool.imap(process_single_doi, args_list),
            total=len(args_list),
            desc="处理DOI"
        ))
        
    # 处理结果
    for doi, success, error_msg in results:
        if success:
            success_dois.append(doi)
        else:
            failed_dois.append((doi, error_msg))
            print(f"处理失败 {doi}: {error_msg}")
    
    # 保存处理结果
    print(f"\n处理完成:")
    print(f"成功: {len(success_dois)} 个DOI")
    print(f"失败: {len(failed_dois)} 个DOI")
    
    # 保存失败的DOI和错误信息
    if failed_dois:
        failed_df = pd.DataFrame(failed_dois, columns=['DOI', 'Error'])
        failed_df.to_csv('failed_dois.csv', index=False)
        print("失败的DOI已保存到 failed_dois.csv")

if __name__ == '__main__':
    main()
