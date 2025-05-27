import time
import requests
import base64
import os
import os
os.environ['DISPLAY'] = ':0'
os.environ['XAUTHORITY'] = '/home/tju/.Xauthority'
import pyautogui
from typing import List, Dict, Any
from dataclasses import dataclass
from typing import Optional
from enum import Enum
import datetime

def encode_image(image_path):
    """将图片转换为base64格式"""
    with open(image_path, "rb") as f:
        return base64.b64encode(f.read()).decode('utf-8')

def save_base64_image(base64_str, output_path):
    """将base64图片保存为文件"""
    img_data = base64.b64decode(base64_str)
    with open(output_path, 'wb') as f:
        f.write(img_data)

def check_pdf_num(path='~/Downloads'):
    path = os.path.expanduser(path)
    return len(os.listdir(path))

def get_newest_file(path='~/Downloads'):
    path = os.path.expanduser(path)
    files = [f for f in os.listdir(path)]
    if not files:
        return None
    newest_file = max(files, key=lambda x: os.path.getctime(os.path.join(path, x)))
    return newest_file


def parse_image(image_path, omniparser_url="http://localhost:8000",box_threshold=0.01,
    iou_threshold=0.1,
    use_paddleocr=True,):
    """调用omniparser服务解析图片"""
    try:
        # 检查服务是否在运行
        probe_response = requests.get(f"{omniparser_url}/probe/")
        if probe_response.status_code != 200:
            raise Exception("Omniparser service is not available")

        # 转换图片为base64
        image_base64 = encode_image(image_path)
        
        # 调用解析服务
        response = requests.post(
            f"{omniparser_url}/parse/",
            json={"base64_image": image_base64,
                  "params":
                      {
                        "box_threshold": box_threshold,
                        "iou_threshold": iou_threshold,
                        "use_paddleocr": use_paddleocr,
                      }
                      }
        )

        if response.status_code != 200:
            raise Exception(f"Parse failed with status code: {response.status_code}")
            
        return response.json()
    
    except requests.exceptions.RequestException as e:
        print(f"Error connecting to omniparser service: {e}")
        return None
    
    
class DownloadStep(Enum):
    BROWSER_DETECTION = "browser_detection"
    INITIAL_NAVIGATION = "initial_navigation"
    PAGE_LOAD_CHECK = "page_load_check"
    PDF_BUTTON_SEARCH = "pdf_button_search"
    PDF_DOWNLOAD_CHECK = "pdf_download_check"
    VERIFICATION_HANDLER = "verification_handler"
    DOWNLOAD_COMPLETION = "download_completion"
    
@dataclass
class DownloadState:
    doi: str
    initial_pdf_count: int
    current_step: DownloadStep = DownloadStep.BROWSER_DETECTION
    max_steps: int = 10
    step_count: int = 0
    last_screen: Optional[Dict] = None
    current_url: Optional[str] = None
    publisher: Optional[str] = None
    previous_action: Optional[str] = None
    challenge_type: Optional[str] = None
    verification_required: bool = False
    browser_address_bar_coords: Optional[List[float]] = None

    
    
def grab_screen(state: DownloadState) -> Dict[str, Any]:
    """截取屏幕并分析按钮"""
    current_datetime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
    image_path = f"./tmp/screen_{state.current_step.value}_{current_datetime}.png"
    pyautogui.screenshot(image_path)
    
    # 使用OmniParser服务解析图片
    result = parse_image(image_path)
    if not result:
        return {"error": "Failed to parse screen"}
    
    output_path = f"./tmp/labeled_{state.current_step.value}_{current_datetime}.png"
    save_base64_image(result["som_image_base64"], output_path)
    
    # 更新状态
    
    return {
        "image_path": image_path,
        "labeled_image_path": output_path,
        "elements": '\n'.join([f"icon {i}: {content}" for i, content in enumerate(result["parsed_content_list"])]),
        "som_image_base64": result["som_image_base64"]
    }
    

def execute_command(action: str, coordinates: Optional[List[float]] = None, url: str = None):
    """执行鼠标和键盘命令
    coordinates: [x1, y1, x2, y2]
    text: url
    """
    # if coordinates is None and action != "wait":
    #     raise ValueError("Coordinates are required for non-wait actions")
    screen_width, screen_height = pyautogui.size()

    
    if action == "search_pdf_button":
        return "Search PDF button clicked"
    
    if coordinates is not None:
        # 提取归一化坐标
        x1_norm, y1_norm, x2_norm, y2_norm = coordinates

        # 转换为实际像素位置
        x1 = int(x1_norm * screen_width)
        y1 = int(y1_norm * screen_height)
        x2 = int(x2_norm * screen_width)
        y2 = int(y2_norm * screen_height)

        # 计算中心点
        center_x = int((x1 + x2) / 2)
        center_y = int((y1 + y2) / 2)
    
    if action == "move_and_left_click":
        if coordinates is None:
            raise ValueError("Coordinates are required for left click action")
        pyautogui.moveTo(center_x, center_y)
        pyautogui.click()
        # time.sleep(2)
        return f"Clicked at ({center_x}, {center_y})"
    if action == "wait":
        time.sleep(2)
        return "Waited for 2 seconds"
    if action == "move_and_enter_url":
        if coordinates is None:
            raise ValueError("Coordinates are required for enter url action")
        # if doi is None:
        #     raise ValueError("Url is required for enter url action")
        # url = f"https://doi.org/{doi}"
        pyautogui.moveTo(center_x, center_y)
        pyautogui.click()
        pyautogui.hotkey('ctrl', 'a')
        pyautogui.press('backspace')
        pyautogui.typewrite(url, interval=0.03)
        pyautogui.press('enter')
        time.sleep(4)
        return f"Entered URL: {url}"
    if action == "page_down":
        pyautogui.press("pagedown")
        time.sleep(2)
        return "Pressed page down"
    if action == "page_up":
        pyautogui.press("pageup")
        time.sleep(2)
        return "Pressed page up"
    else:
        raise ValueError(f"Unsupported action: {action}")
    
    
        
    