import pandas as pd
import os
from openai import OpenAI
import base64
import json
from typing import Optional, Dict, Any
import time

def encode_image(image_path: str) -> str:
    """Encode image to base64 string."""
    if not os.path.exists(image_path):
        raise FileNotFoundError(f"Image file not found: {image_path}")
    try:
        with open(image_path, "rb") as image_file:
            return base64.b64encode(image_file.read()).decode("utf-8")
    except Exception as e:
        raise Exception(f"Failed to encode image: {str(e)}")

class LLMCaller:
    def __init__(self, 
                 api_key: Optional[str] = None,
                 model: str = "deepseek-chat",
                 base_url: str = "api.deepseek.com/v1",
                 max_retries: int = 3,
                 retry_delay: float = 2.0):
        """Initialize LLM caller with configuration."""
        if api_key:
            os.environ["OPENAI_API_KEY"] = api_key
        if base_url:
            self.client = OpenAI(base_url=base_url)
        else:
            self.client = OpenAI()
        self.model = model
        self.max_retries = max_retries
        self.retry_delay = retry_delay

    def call_llm(self, 
                prompt: str, 
                system_message: str = "You are a helpful assistant.",
                response_json: bool = False,
                stream: bool = False) -> str:
        """Call LLM with prompt and return response with retry mechanism."""
        last_error = None
        for attempt in range(self.max_retries):
            try:
                kwargs = {
                    "model": self.model,
                    "messages": [
                        {"role": "system", "content": system_message},
                        {"role": "user", "content": prompt},
                    ],
                    "stream": stream,
                    "temperature": 1.0  # 添加temperature参数，0.2适合结构化输出

                }
                
                if response_json:
                    kwargs["response_format"] = {"type": "json_object"}
                
                completion = self.client.chat.completions.create(**kwargs)
                
                if stream:
                    response_text = ""
                    for chunk in completion:
                        if chunk.choices and chunk.choices[0].delta.content is not None:
                            chunk_content = chunk.choices[0].delta.content
                            response_text += chunk_content
                            # print(chunk_content, end="", flush=True)
                    print()  # New line after streaming completes
                    return response_text
                else:
                    return completion.choices[0].message.content
                    
            except Exception as e:
                last_error = e
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay)
                    continue
        
        raise Exception(f"LLM API call failed after {self.max_retries} attempts. Last error: {str(last_error)}")

class VisualLLMCaller:
    def __init__(self, 
                 model: str = "qwen2.5-vl-72b-instruct",
                 api_key: Optional[str] = None,
                 base_url: Optional[str] = None,
                 max_retries: int = 3,
                 retry_delay: float = 2.0):
        """Initialize Visual LLM caller with configuration."""
        if api_key:
            os.environ["OPENAI_API_KEY"] = api_key
        if base_url:
            self.client = OpenAI(base_url=base_url)
        else:
            self.client = OpenAI()
        self.model = model
        self.max_retries = max_retries
        self.retry_delay = retry_delay

    def call_llm(self, 
                 prompt: str, 
                 images_list: str,
                 system_message: str = "You are a helpful assistant.",
                 response_json: bool = False,
                 stream: bool = False) -> str:
        """Call Visual LLM with prompt and image, return response with retry mechanism."""
        # if not os.path.exists(images_path):
        #     raise FileNotFoundError(f"Image file not found: {images_path}")
        
        # 构建消息内容
        content = []
        images_name = []
        
        # 首先添加文本提示
        content.append({
            "type": "text",
            "text": prompt
        })
        
        # 添加图片
        for image_file in images_list:
            if image_file.endswith(('.png', '.jpg', '.jpeg')):
                content.append({
                    "type": "image_url",
                    "image_url": {
                        "url": f"data:image/png;base64,{encode_image(os.path.join(images_path, image_file))}"
                    }
                })
                images_name.append(image_file)
            
        last_error = None
        for attempt in range(self.max_retries):
            try:
                kwargs = {
                    "model": self.model,
                    "messages": [{"role": "user", "content": content}],
                    "stream": stream,
                    "top_p": 0.20
                }

                if response_json:
                    kwargs["response_format"] = {"type": "json_object"}

                completion = self.client.chat.completions.create(**kwargs)
                
                if stream:
                    response_text = ""
                    for chunk in completion:
                        if chunk.choices and chunk.choices[0].delta.content is not None:
                            chunk_content = chunk.choices[0].delta.content
                            response_text += chunk_content
                    print()
                    return response_text
                else:
                    return completion.choices[0].message.content
                    
            except Exception as e:
                last_error = e
                if attempt < self.max_retries - 1:
                    time.sleep(self.retry_delay)
                    continue
        
        raise Exception(f"Visual LLM API call failed after {self.max_retries} attempts. Last error: {str(last_error)}")
