import time
from openai import OpenAI
import os
from llm_tools import get_newest_file, grab_screen, execute_command, DownloadState, check_pdf_num
import json
from dataclasses import dataclass
from typing import Any, List, Dict, Optional
from enum import Enum
from prompts import DownloadPrompts
import re

class DownloadStep(Enum):
    INITIAL_NAVIGATION = "initial_navigation"
    PAGE_LOAD_CHECK = "page_load_check"
    PDF_BUTTON_SEARCH = "pdf_button_search"
    PDF_DOWNLOAD_CHECK = "pdf_download_check"
    VERIFICATION_HANDLER = "verification_handler"
    DOWNLOAD_COMPLETION = "download_completion"

# 设置X显示
os.environ['DISPLAY'] = ':0'
os.environ['XAUTHORITY'] = '/home/tju/.Xauthority'


# LLM_CONFIG = {
#     "model_name": os.getenv("MODEL_NAME", "gpt-4o"),
#     "base_url": os.getenv("BASE_URL", "https://ai98.vip/v1"),
#     "api_key": os.getenv("API_KEY", "sk-ATTA33l8XGlDgasiKjXKFEm2EF0lm0hreKesGB6e8dqoZxct")
# }

LLM_CONFIG = {
    "model_name": os.getenv("MODEL_NAME", "qwen-vl-max"),
    "base_url": os.getenv("BASE_URL", "https://dashscope.aliyuncs.com/compatible-mode/v1"),
    "api_key": os.getenv("API_KEY", "sk-7fbb65e9ab6f40bcb3acc916288a324b")
}

@dataclass
class DownloadState:
    doi: str
    initial_pdf_count: int
    current_step: DownloadStep = DownloadStep.INITIAL_NAVIGATION
    max_steps: int = 10
    step_count: int = 0
    current_url: Optional[str] = None
    publisher: Optional[str] = None
    previous_action: Optional[str] = None
    challenge_type: Optional[str] = None
    verification_required: bool = False
    browser_address_bar_coords: Optional[List[float]] = None

def parse_response(response_text: str) -> Dict[str, Any]:
    """解析LLM的响应"""
    try:
        # 尝试直接解析
        result = json.loads(response_text)
        return result
    except json.JSONDecodeError:
        # 如果直接解析失败，尝试从文本中提取JSON
        try:
            # 查找JSON代码块
            json_match = re.search(r'```json\s*(.*?)\s*```', response_text, re.DOTALL)
            if json_match:
                json_str = json_match.group(1).strip()
                # 替换单引号为双引号（处理常见错误）
                json_str = re.sub(r"'([^']*)':", r'"\1":', json_str)  # 替换键中的单引号
                json_str = re.sub(r':\s*\'([^\']*)\'', r': "\1"', json_str)  # 替换值中的单引号
                return json.loads(json_str)
            else:
                # 如果没有找到JSON代码块，尝试查找第一个{和最后一个}之间的内容
                json_match = re.search(r'{.*}', response_text, re.DOTALL)
                if json_match:
                    json_str = json_match.group(0)
                    # 替换单引号为双引号
                    json_str = re.sub(r"'([^']*)':", r'"\1":', json_str)
                    json_str = re.sub(r':\s*\'([^\']*)\'', r': "\1"', json_str)
                    return json.loads(json_str)
                else:
                    raise ValueError("无法从响应中提取JSON")
        except Exception as e:
            print(f"解析JSON失败: {e}")
            print(f"原始响应: {response_text}")
            # 尝试手动解析
            try:
                # 提取关键信息
                is_browser_match = re.search(r'"is_browser":\s*(true|false)', response_text)
                action_match = re.search(r'"action":\s*[\'"]([^\'"]*)[\'"]', response_text)
                bbox_match = re.search(r'"icon_bbox":\s*\[([\d\., ]+)\]', response_text)
                
                result = {}
                if is_browser_match:
                    result["is_browser"] = is_browser_match.group(1) == "true"
                if action_match:
                    result["action"] = action_match.group(1)
                if bbox_match:
                    bbox_str = bbox_match.group(1)
                    result["icon_bbox"] = [float(x.strip()) for x in bbox_str.split(',')]
                
                if result:
                    print("手动解析结果:", result)
                    return result
            except Exception as e2:
                print(f"手动解析也失败: {e2}")
            
            # 如果所有解析方法都失败，返回默认值
            return {"action": "wait", "icon_bbox": None, "url": None}

class DownloadAgent:
    def __init__(self, doi: str):
        self.state = DownloadState(
            doi=doi,
            initial_pdf_count=check_pdf_num()
        )
        self.llm = self._init_llm()
        self.prompts = DownloadPrompts()
        self.current_screen = grab_screen(self.state)
        self.current_prompt = self._get_current_prompt()
        self.logs = []
    
    def log(self, msg):
        self.logs.append(msg)

    def _get_current_prompt(self) -> List[Dict[str, Any]]:
        """获取当前步骤的提示，包含屏幕元素信息"""
        prompt_map = {
            DownloadStep.INITIAL_NAVIGATION: self.prompts.INITIAL_NAVIGATION,
            DownloadStep.PAGE_LOAD_CHECK: self.prompts.PAGE_LOAD_CHECK,
            DownloadStep.PDF_BUTTON_SEARCH: self.prompts.PDF_BUTTON_SEARCH,
            DownloadStep.PDF_DOWNLOAD_CHECK: self.prompts.PDF_DOWNLOAD_CHECK,
            DownloadStep.VERIFICATION_HANDLER: self.prompts.VERIFICATION_HANDLER,
            DownloadStep.DOWNLOAD_COMPLETION: self.prompts.DOWNLOAD_COMPLETION
        }
        prompt_template = prompt_map[self.state.current_step]
        text_message ={"type": "text",
                       "text":   prompt_template.format(
                           doi=self.state.doi,
                           current_url=self.state.current_url,
                           publisher=self.state.publisher,
                           current_step=self.state.current_step.value,
                           previous_action=self.state.previous_action,
                           challenge_type=self.state.challenge_type,
                           initial_count=self.state.initial_pdf_count,
                           current_count=check_pdf_num(),
                           address_bar_coords=self.state.browser_address_bar_coords,
                           elements_info=self.current_screen.get("elements", []),
                           som_image_base64=self.current_screen.get("som_image_base64", "")
                       )}
        image_message = {"type": "image_url", 
                         'image_url': {'url': f"data:image/png;base64,{self.current_screen.get('som_image_base64', '')}"}}
        message = [{"role": "user", "content": [text_message, image_message]}]
        return message

    def _init_llm(self):
        """初始化LLM"""
        chat_llm = OpenAI(
            api_key=LLM_CONFIG["api_key"],
            base_url=LLM_CONFIG["base_url"],
        )
        return chat_llm
    

    def update_state(self, result: Dict[str, Any]):
        """更新状态"""
        self.state.last_screen = self.current_screen
        self.current_prompt = self._get_current_prompt()
        if result.get("is_browser") == True:
            self.state.current_step = DownloadStep.PAGE_LOAD_CHECK
            self.current_prompt = self._get_current_prompt()
        if result.get("page_loaded") == True:
            self.state.current_step = DownloadStep.PDF_BUTTON_SEARCH
            self.current_prompt = self._get_current_prompt()
        if result.get("pdf_button_found") == True:
            self.state.current_step = DownloadStep.PDF_DOWNLOAD_CHECK
            self.current_prompt = self._get_current_prompt()

    

    
    def run(self):
        # 1. 先判断当前的截图是否是浏览器
            # 询问
        #循环执行，直到self.state.current_step == DownloadStep.PDF_BUTTON_SEARCH
        self.log("Current step:"+self.state.current_step.value)
        while self.state.current_step != DownloadStep.PAGE_LOAD_CHECK:
            
            self.state.step_count = 0
        
            response = self.llm.chat.completions.create(
                model=LLM_CONFIG["model_name"],
                messages=self.current_prompt,
                stream=True,
                # response_format={"type": "json_object"},
            )
            
            # print(self.current_prompt)
            
            response_text = ""
            for chunk in response:
                if chunk.choices and chunk.choices[0].delta.content is not None:
                    response_text += chunk.choices[0].delta.content
                    # print(chunk.choices[0].delta.content, end="")
            
            self.log(f"data:image/png;base64,"+self.current_screen.get("som_image_base64", "")+"[IMG_END]"+response_text)
            
            # save txt
            # with open("response.txt", "a") as f:
            #     f.write(f"data:image/png;base64,"+self.current_screen.get("som_image_base64", "")+response_text)
            
            # print(self.current_screen.get("image_path"))
            
            # 执行
            result = parse_response(response_text)
            execute_command(result["action"], result["icon_bbox"], f"https://doi.org/{self.state.doi}")
            time.sleep(2)
            
            self.current_screen = grab_screen(self.state)
            
            self.update_state(result)
            
            self.state.step_count += 1
            if self.state.step_count >= self.state.max_steps:
                raise Exception("Max steps reached, aborting...")
        
        # 进行下一个步骤的循环
        self.log("Current step:"+self.state.current_step.value)
        while self.state.current_step != DownloadStep.PDF_BUTTON_SEARCH:
            
            self.state.step_count = 0
            
            response = self.llm.chat.completions.create(
                model=LLM_CONFIG["model_name"],
                messages=self.current_prompt,
                stream=True,
                # response_format={"type": "json_object"},
            )
            
            response_text = ""
            for chunk in response:
                if chunk.choices and chunk.choices[0].delta.content is not None:
                    response_text += chunk.choices[0].delta.content
                    # print(chunk.choices[0].delta.content, end="")
            
            self.log(f"data:image/png;base64,"+self.current_screen.get("som_image_base64", "")+"[IMG_END]"+response_text)

            result = parse_response(response_text)
            execute_command(result["action"], result["icon_bbox"], f"https://doi.org/{self.state.doi}")
            time.sleep(2)

            self.current_screen = grab_screen(self.state)
            
            self.update_state(result)
            
            self.state.step_count += 1
            if self.state.step_count >= self.state.max_steps:
                raise Exception("Max steps reached, aborting...")
        
        while self.state.current_step != DownloadStep.PDF_DOWNLOAD_CHECK:
            self.log("Current step:"+self.state.current_step.value)
            
            self.state.step_count = 0
            response = self.llm.chat.completions.create(
                model=LLM_CONFIG["model_name"],
                messages=self.current_prompt,
                stream=True,
                # response_format={"type": "json_object"},
            )
            response_text = ""
            for chunk in response:
                if chunk.choices and chunk.choices[0].delta.content is not None:
                    response_text += chunk.choices[0].delta.content
                    # print(chunk.choices[0].delta.content, end="")
            
            self.log(f"data:image/png;base64,"+self.current_screen.get("som_image_base64", "")+"[IMG_END]"+response_text)

            result = parse_response(response_text)
            
            if result.get("pdf_button_found") == False:
                raise Exception("PDF button not found, aborting...")
            
            execute_command(result["action"], result["icon_bbox"], f"https://doi.org/{self.state.doi}")
            time.sleep(2)

            self.current_screen = grab_screen(self.state)
            
            self.update_state(result)
            
            self.state.step_count += 1
            if self.state.step_count >= self.state.max_steps:
                raise Exception("Max steps reached, aborting...")
            
        while self.state.current_step != DownloadStep.DOWNLOAD_COMPLETION:
            self.log("Current step:"+self.state.current_step.value)
            self.state.step_count = 0
            
            # 检查PDF数量
            current_count = check_pdf_num()
            if current_count > self.state.initial_pdf_count:
                # 如果最新的
                while True:
                    newst_file = get_newest_file()
                    # 如果newst_file以pdf结尾
                    if newst_file.endswith(".pdf"):
                        self.state.current_step = DownloadStep.DOWNLOAD_COMPLETION
                        break
                    elif newst_file.endswith(".crdownload"):
                        # 等待
                        time.sleep(2)
                        continue
                    else:
                        raise Exception("Unknown file type, aborting...")
            else:
                while self.state.current_step != DownloadStep.DOWNLOAD_COMPLETION:
                    response = self.llm.chat.completions.create(
                        model=LLM_CONFIG["model_name"],
                        messages=self.current_prompt,
                        stream=True,
                        # response_format={"type": "json_object"},
                    )
                    
                    response_text = ""
                    for chunk in response:
                        if chunk.choices and chunk.choices[0].delta.content is not None:
                            response_text += chunk.choices[0].delta.content
                            # print(chunk.choices[0].delta.content, end="")
                    
                    self.log(f"data:image/png;base64,"+self.current_screen.get("som_image_base64", "")+"[IMG_END]"+response_text)

                    result = parse_response(response_text)
                    execute_command(result["action"], result["icon_bbox"], f"https://doi.org/{self.state.doi}")
                    time.sleep(2)
                    self.current_screen = grab_screen(self.state)
                    
                    self.update_state(result)
                    
                    self.state.step_count += 1
                    if self.state.step_count >= self.state.max_steps:
                        self.log("Max steps reached, aborting...")
                        raise Exception("Max steps reached, aborting...")
                    
                    if check_pdf_num() > self.state.initial_pdf_count and "Chrome" not in get_newest_file():
                        
                        self.state.current_step = DownloadStep.DOWNLOAD_COMPLETION
                        break
                    
        self.log("Wait for download completion...")
        while True:
            if check_pdf_num() > self.state.initial_pdf_count and get_newest_file().endswith(".pdf"):
                self.log("Download completed!, PDF file: "+"~/Downloads/" +get_newest_file())
                # os.system(f"mv ~/Downloads/{get_newest_file()} ~/Downloads/{self.state.doi}.pdf")
                break
            time.sleep(2)
            
        
    

if __name__ == "__main__":
    # agent = DownloadAgent("10.1016/j.cej.2025.162723")
    agent = DownloadAgent("10.1002/adma.202405163")
    agent.run()
          
        
        
        
    