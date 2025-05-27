# app.py  –  full version with live console for DownloadAgent
import uuid
from openai import OpenAI
import streamlit as st
from PIL import Image
import os
import threading, io, contextlib, time            # NEW
from pathlib import Path
import streamlit.components.v1 as components
import pyautogui
import datetime
from download_agent import DownloadAgent          # ← your existing class
from test import PrintCollector
from streamlit.runtime.scriptrunner import add_script_run_ctx  # Streamlit ≥1.27
import base64
import tempfile

# ---- no-VNC ---------------------------------------------------------------
NOVNC_URL = (
    "http://172.26.152.0:6080/vnc_auto.html"
    "?password=123@tju"
    "&view_only=1"
    "&autoconnect=true"
    "&resize=scale"
)

# ---- Streamlit page -------------------------------------------------------
st.set_page_config(
    page_title="AI Research Assistant",
    page_icon="🤖",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={

        'About': "# OLED Research Assistant v1.0.0"
    }
)

# ---- global CSS (kept from original) -------------------------------------
st.markdown(
    """
    <style>
    :root {
        --background-color: #ffffff;
        --text-color: #31333F;
        --sidebar-bg: #f8f9fa;
        --card-bg: #ffffff;
        --hover-bg: #e9ecef;
        --border-color: #dee2e6;
        --primary-color: #1f77b4;
    }
    /* Keep necessary styles */
    .agent-title{font-size:20px;font-weight:bold;margin-top:10px;color:var(--text-color);}
    [data-testid="stSidebar"]{background:var(--sidebar-bg);padding:2rem 1rem;}
    .stApp{background:var(--background-color);}
    .stButton>button{background:var(--primary-color);color:#fff;}
    .stCode{max-height:300px;overflow-y:auto;}        /* log panel height */
    /* Optimize image display */
    img {image-rendering: -webkit-optimize-contrast; image-rendering: crisp-edges;}
    [data-testid="caption"] {margin-top: 0.5rem; font-size: 0.9rem; color: #555;}
    .agent-container {
        border: 1px solid #dee2e6;
        border-radius: 10px;
        padding: 20px;
        box-shadow: 0 4px 6px rgba(0,0,0,.1);
        transition: .3s;
        text-align: center;
        background: white;
        height: 100%;
        display: flex;
        flex-direction: column;
        justify-content: space-between;
        min-height: 300px;
        align-items: center;
        position: relative;
    }
    .agent-container:hover {
        transform: scale(1.02);
        box-shadow: 0 6px 8px rgba(0,0,0,.15);
    }
    .agent-icon {
        font-size: 100px;
        margin: 20px 0;
        height: 120px;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    .agent-icon img {
        max-height: 120px;
        max-width: 100%;
        object-fit: contain;
    }
    .agent-title {
        font-size: 20px;
        font-weight: bold;
        margin-bottom: 15px;
        color: #31333F;
    }
    .agent-description {
        flex-grow: 1;
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
        padding: 0 10px;
    }
    .agent-button {
        display: block;
        background: #1f77b4;
        color: white;
        padding: 12px 24px;
        border-radius: 30px;
        text-decoration: none;
        font-weight: bold;
        transition: all 0.3s;
        width: 80%;
        box-sizing: border-box;
        box-shadow: 0 3px 5px rgba(0,0,0,0.2);
        border: none;
        letter-spacing: 1px;
        position: absolute;
        bottom: 20px;
        left: 50%;
        transform: translateX(-50%);
        text-align: center;
    }
    .agent-button:hover {
        background: #1565C0;
        transform: translateY(-2px);
        box-shadow: 0 5px 8px rgba(0,0,0,0.3);
    }
    .agent-button:active {
        transform: translateY(0);
        box-shadow: 0 2px 3px rgba(0,0,0,0.2);
    }
    [data-testid="column"] {
        padding: 0 10px;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# ---- sidebar --------------------------------------------------------------
def sidebar_nav():
    with st.sidebar:
        st.image("assets/logo.png", width=400, channels="BGR")
        st.markdown('<div style="text-align:center;"><h2>AI Research Assistant</h2></div>',
                    unsafe_allow_html=True)

        menu = st.radio("Navigation Menu", ["Home",
                                            "Paper Collection Agent",
                                            "Information Extraction Agent",
                                            "Deep Research Agent"],
                        label_visibility="collapsed")
        st.markdown("---")
        with st.expander("Settings"):
            st.checkbox("Dark Mode")
            st.slider("Font Size", 12, 24, 16)
        st.markdown("---")
        st.write("Version 1.0.0  ")
    return menu

# ---- pages ----------------------------------------------------------------
def show_home():
    st.title("Welcome to OLED Research Assistant")
    st.markdown("### Choose your AI Agent")
    
    col1, col2, col3 = st.columns(3)

    # Use st.image and st.container to replace HTML cards, apply use_container_width =True
    with col1:
        with st.container():
            st.markdown("""
            <div class="agent-container" onclick="parent.streamlitPyToJsMessage({
                    type: 'streamlit:setComponentValue',
                    value: 'Paper Collection Agent',
                    dataType: 'text',
                    key: 'current_page'
                }); window.location.reload();" style="cursor: pointer;">
                <div>
                    <div class="agent-title">Paper Collection Agent</div>
                    <div class="agent-icon">
                        📥
                    </div>
                </div>
                <div class="agent-description">
                    <p>Automatically download research papers from various sources.</p>
                </div>
            </div>
            """, unsafe_allow_html=True)

    with col2:
        with st.container():
            st.markdown("""
            <div class="agent-container" onclick="parent.streamlitPyToJsMessage({
                    type: 'streamlit:setComponentValue',
                    value: 'Information Extraction Agent',
                    dataType: 'text',
                    key: 'current_page'
                }); window.location.reload();" style="cursor: pointer;">
                <div>
                    <div class="agent-title">Information Extraction Agent</div>
                    <div class="agent-icon">
                        📊
                    </div>
                </div>
                <div class="agent-description">
                    <p>Extract key information and summaries from papers.</p>
                </div>
            </div>
            """, unsafe_allow_html=True)

    with col3:
        with st.container():
            st.markdown("""
            <div class="agent-container" onclick="parent.streamlitPyToJsMessage({
                    type: 'streamlit:setComponentValue',
                    value: 'Deep Research Agent',
                    dataType: 'text',
                    key: 'current_page'
                }); window.location.reload();" style="cursor: pointer;">
                <div>
                    <div class="agent-title">Deep Research Agent</div>
                    <div class="agent-icon">
                        🤖
                    </div>
                </div>
                <div class="agent-description">
                    <p>Train models and get AI-powered insights.</p>
                </div>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### Overview")
    col1, col2, col3 = st.columns([1, 3, 1])
    with col2:
        # Use PNG version of the image, potentially better quality
        st.image("assets/Figure1.png", caption='Overview of the AI Agents system', use_container_width=True)
        
        # If the above method doesn't work, try uncommenting one of the methods below:
        # st.image("assets/Figure1.jpg", caption='Overview of the AI Agents system', use_container_width =True)
        # st.image("assets/Figure1.jpg", caption='Overview of the AI Agents system', width=800)
        
        # Or use HTML to directly insert the image, can better control the display of the image
        # st.markdown("""
        #     <div style="text-align:center;">
        #         <img src="assets/Figure1.png" style="max-width:100%; height:auto; image-rendering:high-quality;" />
        #         <p style="margin-top:5px; font-style:italic; color:#555;">Overview of the AI Agents system</p>
        #     </div>
        # """, unsafe_allow_html=True)

def process_logs_to_html(logs, with_avatar=True):
    """
    Process logs into HTML format, supporting text and image display, with optional avatars
    
    Args:
        logs: List containing text and base64 encoded images
        with_avatar: Whether to add avatars
        
    Returns:
        str: Processed HTML content
    """
    # Define avatar style and HTML
    avatar_html = ""
    if with_avatar:
        avatar_html = """
        <div style="display: flex; align-items: flex-start; margin-bottom: 10px;">
            <div style="width: 36px; height: 36px; border-radius: 50%; background-color: #ffffff; 
                        display: flex; justify-content: center; align-items: center; margin-right: 10px;
                        flex-shrink: 0; color: #1f77b4; font-weight: bold; border: 1px solid #dee2e6;">
                🤖
            </div>
            <div style="flex-grow: 1; background-color: #ffffff; border-radius: 10px; padding: 10px; 
                      border: 1px solid #dee2e6; box-shadow: 0 1px 3px rgba(0,0,0,0.1);">
        """
    
    # Start building content
    html_content = avatar_html if with_avatar else ""
    html_content += "<pre style='margin: 0; white-space: pre-wrap;'>"
    
    for log in logs:
        # Check if contains base64 image
        if "data:image/png;base64," in log and "[IMG_END]" in log:
            # Use separators to precisely split
            parts = log.split("data:image/png;base64,", 1)
            prefix_text = parts[0]  # Text before image
            
            # Split base64 data and subsequent text
            img_and_text = parts[1].split("[IMG_END]", 1)
            base64_data = img_and_text[0]  # base64 image data
            
            # Add prefix text
            html_content += prefix_text + "</pre>"
            
            # Add image
            html_content += f'<img src="data:image/png;base64,{base64_data}" style="max-width:50%; margin:10px 0; border-radius: 5px;">'
            
            # Add suffix text (if any)
            html_content += "<pre style='margin: 0; white-space: pre-wrap;'>"
            if len(img_and_text) > 1:
                text_after = img_and_text[1]
                html_content += text_after
        else:
            # Regular text log
            html_content += log + "\n"
        
    html_content += "</pre>"
    
    # If avatar was added, close div
    if with_avatar:
        html_content += "</div></div>"
    
    return html_content

def show_paper_download():
    st.title("Paper Collection Agent")
    st.markdown("✨ This agent is used to download scientific papers using browsers. A remote desktop is provided for you to operate.")
    
    # Add video autoplay - Method 2: Using HTML (supports autoplay)
    video_path = "assets/video-1080.webm"
    
    # Create a unique ID to reference this video element on the page
    video_id = "autoplay-video"
    
    # Use HTML and JavaScript to implement autoplay
    st.markdown(f"""
    <div style="display: flex; justify-content: center; margin: 20px 0;">
        <video id="{video_id}" autoplay loop muted playsinline 
               style="max-width: 100%; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.2);">
            <source src="data:video/webm;base64,{base64.b64encode(open(video_path, 'rb').read()).decode()}" type="video/webm">
            Your browser does not support the video tag.
        </video>
    </div>
    <script>
        // Ensure video autoplays
        document.addEventListener('DOMContentLoaded', (event) => {{
            document.getElementById('{video_id}').play();
        }});
    </script>
    """, unsafe_allow_html=True)
    
    st.markdown("---")

    col1, col2 = st.columns([3, 2])
    with col1:
    # ---- noVNC remote desktop --------------------------------------------
        st.subheader("Remote Desktop (noVNC)")
        components.iframe(NOVNC_URL, height=800, scrolling=True)

    with col2:
    # ---------------- Search Area -----------------
        doi = st.text_input("Enter DOI", value="")
        # output_placeholder = st.empty().container(height=500)
        output_placeholder = st.container(height=700)
        

        # Initial welcome message
        welcome_html = process_logs_to_html(["👋 Hi, I am an AI agent for downloading scientific papers using browsers. Please enter the DOI and click 'Download' to start."], with_avatar=True)
        output_placeholder.markdown(welcome_html, unsafe_allow_html=True)
        
        if st.button("Download") and doi:
            agent = DownloadAgent(doi)
            
            def background_task():
                agent.run()
            thread = threading.Thread(target=background_task)
            thread.start()
            
            # Main thread refreshes display content
            logs_length = 0
            while thread.is_alive():
                if len(agent.logs) > logs_length:
                    new_logs = agent.logs[logs_length:]  # Get all new logs
                    logs_length = len(agent.logs)  # Update log length
                    
                    # Process log content and display
                    html_content = process_logs_to_html(new_logs, with_avatar=True)
                    output_placeholder.markdown(html_content, unsafe_allow_html=True)
                    
                time.sleep(0.2)

            # Check one more time
            if len(agent.logs) > logs_length:
                new_logs = agent.logs[logs_length:]  # Get all new logs
                logs_length = len(agent.logs)  # Update log length
                
                # Process log content and display
                html_content = process_logs_to_html(new_logs, with_avatar=True)
                output_placeholder.markdown(html_content, unsafe_allow_html=True)
                
            # Add completion message
            completion_html = process_logs_to_html(["✅ Download completed."], with_avatar=True)
            output_placeholder.markdown(completion_html, unsafe_allow_html=True)

# ---------------------------------------------------------------------------
def show_paper_parser():
    st.title("Information Extraction Agent")
    
    # Add Token input
    if 'mineru_token' not in st.session_state:
        st.session_state.mineru_token = ""
    
    if 'deepseek_api_key' not in st.session_state:
        st.session_state.deepseek_api_key = ""
    
    with st.expander("API Settings", expanded=False):
        col1, col2 = st.columns(2)
        
        with col1:
            token_input = st.text_input(
                "MinerU API Token (https://mineru.net/apiManage/token)", 
                value=st.session_state.mineru_token,
                type="password",
                help="Please enter your MinerU API token for PDF parsing"
            )
            if token_input != st.session_state.mineru_token:
                st.session_state.mineru_token = token_input
                st.success("MinerU Token updated")
        
        with col2:
            deepseek_token = st.text_input(
                "DeepSeek API Key (https://deepseek.com/api)", 
                value=st.session_state.deepseek_api_key,
                type="password",
                help="Please enter your DeepSeek API key for structured data extraction"
            )
            if deepseek_token != st.session_state.deepseek_api_key:
                st.session_state.deepseek_api_key = deepseek_token
                st.success("DeepSeek API Key updated")
    
    # File upload area
    uploaded_file = st.file_uploader("Upload PDF File", type="pdf")
    
    if uploaded_file:
        # Display file information
        file_details = {
            "Filename": uploaded_file.name,
            "File Type": uploaded_file.type,
            "File Size": f"{uploaded_file.size / 1024:.2f} KB"
        }
        st.json(file_details)
        
        # Create temporary folder to store parsing results
        if 'temp_dir' not in st.session_state:
            st.session_state.temp_dir = tempfile.mkdtemp()
            
        # Parsing settings
        with st.expander("Parsing Settings", expanded=True):
            extract_smiles = st.checkbox("Extract molecule SMILES (may take longer)", value=False)
            extract_tables = st.checkbox("Extract table data", value=True)
            extract_formulas = st.checkbox("Extract mathematical formulas", value=True)
            use_deepseek = st.checkbox("Use DeepSeek AI for structured data extraction", value=True)
            
        # Parse button
        if st.button("Start Parsing", key="extract_button"):
            if not st.session_state.mineru_token:
                st.error("Please configure MinerU API Token in settings first")
            elif use_deepseek and not st.session_state.deepseek_api_key:
                st.error("Please configure DeepSeek API Key in settings to use DeepSeek AI")
            else:
                # Save uploaded file to temp directory
                with tempfile.NamedTemporaryFile(delete=False, suffix='.pdf') as tmp_file:
                    tmp_file.write(uploaded_file.getbuffer())
                    pdf_path = tmp_file.name
                
                try:
                    # Call parsing function
                    from parser_agent import parser_pdf, extract_info
                    results = parser_pdf(
                        pdf_path, 
                        token=st.session_state.mineru_token,
                        output_dir=st.session_state.temp_dir
                    )
                    
                    if results:
                        # Extract structured data using DeepSeek if enabled
                        if use_deepseek and 'markdown_content' in results:
                            st.info("Extracting structured data using DeepSeek AI...")
                            
                            # Prepare output path for DeepSeek results
                            deepseek_output_path = os.path.join(
                                results['output_dir'], 
                                f"{results['data_id']}_deepseek.json"
                            )
                            
                            # Call extract_info function
                            deepseek_results = extract_info(
                                markdown_content=results['markdown_content'],
                                api_key=st.session_state.deepseek_api_key,
                                output_path=deepseek_output_path
                            )
                            
                            if deepseek_results:
                                results['deepseek_results'] = deepseek_results
                                results['deepseek_path'] = deepseek_output_path
                        
                        st.session_state.parsing_results = results
                        
                        # Display parsing results
                        st.subheader("Parsing Results")
                        
                        # Create tabs to display different types of content
                        tab_titles = ["Markdown Content", "Images"]
                        if 'deepseek_results' in results:
                            tab_titles.insert(0, "OLED Structured Data")
                            
                        tabs = st.tabs(tab_titles)
                        
                        tab_index = 0
                        if 'deepseek_results' in results:
                            with tabs[0]:  # OLED Structured Data (formerly DeepSeek Analysis)
                                st.subheader("OLED Structured Data")
                                st.json(results['deepseek_results'])
                                
                                # Add download button for DeepSeek results
                                with open(results['deepseek_path'], 'r', encoding='utf-8') as f:
                                    deepseek_content = f.read()
                                    st.download_button(
                                        label="Download Structured Data",
                                        data=deepseek_content,
                                        file_name=f"{results['data_id']}_structured.json",
                                        mime="application/json"
                                    )
                            tab_index += 1
                            
                        with tabs[tab_index]:  # Markdown Content
                            if 'markdown_content' in results:
                                st.markdown(results['markdown_content'])
                            else:
                                st.warning("No Markdown content found")
                            tab_index += 1
                        
                        with tabs[tab_index]:  # Images
                            # Display extracted images
                            if 'images' in results and results['images']:
                                st.write(f"Found {len(results['images'])} images")
                                
                                # Create image gallery with 3 columns
                                cols = st.columns(3)
                                for i, img_info in enumerate(results['images']):
                                    col_idx = i % 3
                                    with cols[col_idx]:
                                        try:
                                            image_path = img_info['path']
                                            image_name = img_info['name']
                                            if os.path.exists(image_path):
                                                st.image(image_path, caption=image_name, use_container_width=True)
                                                with open(image_path, "rb") as img_file:
                                                    img_bytes = img_file.read()
                                                    st.download_button(
                                                        label=f"Download",
                                                        data=img_bytes,
                                                        file_name=image_name,
                                                        mime=f"image/{image_path.split('.')[-1]}"
                                                    )
                                        except Exception as e:
                                            st.error(f"Error displaying image: {str(e)}")
                            else:
                                st.warning("No images found in the document")
                    
                finally:
                    # Clean up temporary PDF file
                    if os.path.exists(pdf_path):
                        os.unlink(pdf_path)

def show_machine_learning():
    st.title("Deep Research Agent")
    
    # 创建标签页
    tab1, tab2 = st.tabs(["OLED Laboratory Preparation Assistant", "Machine Learning Model Training"])
    
    with tab1:
        st.header("🧪 Organic Light-Emitting Diode (OLED) Laboratory Assistant")
        
        # 侧边栏 - API Key 设置
        with st.container():
            st.subheader("Settings")
            col1, col2 = st.columns(2)
            
            with col1:
                # 从环境变量或会话状态获取 API Key
                if "dashscope_api_key" not in st.session_state:
                    st.session_state.dashscope_api_key = os.environ.get("DASHSCOPE_API_KEY", "")
                
                # API Key 输入
                api_key = st.text_input(
                    "DashScope API Key", 
                    value=st.session_state.dashscope_api_key,
                    type="password",
                    help="Enter your DashScope API Key. You can get it from DashScope console: https://dashscope.console.aliyun.com/"
                )
                
                # 保存 API Key 到会话状态
                if api_key != st.session_state.dashscope_api_key:
                    st.session_state.dashscope_api_key = api_key
                    # 重置会话状态
                    if "oled_messages" in st.session_state:
                        st.session_state.oled_messages = []
            
            with col2:
                # 知识库ID输入
                if "knowledge_base_id" not in st.session_state:
                    st.session_state.knowledge_base_id = "lbuni0sw84"
                
                knowledge_base_id = st.text_input(
                    "Knowledge Base ID", 
                    value=st.session_state.knowledge_base_id,
                    help="Enter your knowledge base ID"
                )
                
                if knowledge_base_id != st.session_state.knowledge_base_id:
                    st.session_state.knowledge_base_id = knowledge_base_id
                    # 重置会话状态
                    if "assistant_id" in st.session_state:
                        del st.session_state.assistant_id
                        del st.session_state.thread_id
                        if "oled_messages" in st.session_state:
                            st.session_state.oled_messages = []
            
        col1, col2 = st.columns(2)
        with col1:
            # 模型选择
            model = st.selectbox(
                "Select Model",
                ["qwen-plus", "qwen-max", "qwen-turbo"],
                index=0,
                help="Select the DashScope model to use"
            )
        
        with col2:
            # 分子指纹类型选择
            fp_type = st.selectbox(
                "Molecular Fingerprint Type",
                ["morgan", "maccs"],
                index=0,
                help="Select fingerprint type for molecular similarity search"
            )
            
            # 相似分子返回数量
            top_n = st.slider(
                "Number of Similar Molecules", 
                min_value=1, 
                max_value=20, 
                value=5,
                help="Set the number of similar molecules to return"
            )
        
        # 验证 API Key
        if not st.session_state.dashscope_api_key:
            st.warning("Please enter your DashScope API Key in the settings to continue using the assistant.")
        else:
            # 设置 API Key 到环境变量
            os.environ["DASHSCOPE_API_KEY"] = st.session_state.dashscope_api_key
            
            # 导入需要的库
            try:
                from dashscope import Assistants, Messages, Runs, Threads
                import pandas as pd
                import numpy as np
                import json
                
                # 检查是否已有 rdkit，没有则提示用户安装
                try:
                    from rdkit import Chem
                    from rdkit.Chem import AllChem
                    from rdkit.Chem import MACCSkeys
                    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
                    
                    # 分子指纹相关函数
                    def calculate_morgan_fingerprint(smiles, radius=2, nBits=2048):
                        """Calculate Morgan fingerprint for a molecule"""
                        try:
                            mol = Chem.MolFromSmiles(smiles)
                            if mol is not None:
                                # Use new MorganGenerator instead of deprecated GetMorganFingerprintAsBitVect
                                morgan_gen = GetMorganGenerator(radius=radius, fpSize=nBits)
                                fp = morgan_gen.GetFingerprint(mol)
                                return np.array(fp)
                            else:
                                return np.zeros(nBits)
                        except:
                            return np.zeros(nBits)

                    def calculate_maccs_fingerprint(smiles):
                        """Calculate MACCS fingerprint for a molecule"""
                        try:
                            mol = Chem.MolFromSmiles(smiles)
                            if mol is not None:
                                fp = MACCSkeys.GenMACCSKeys(mol)
                                return np.array(fp)
                            else:
                                return np.zeros(167)  # MACCS fingerprint length is 167
                        except:
                            return np.zeros(167)

                    def search_similar_molecules(query_smiles, fp_type='morgan', top_n=5):
                        """Search for similar molecules based on fingerprint similarity"""
                        # Calculate fingerprint for query molecule
                        if fp_type == 'morgan':
                            query_fp = calculate_morgan_fingerprint(query_smiles)
                            fp_col = 'morgan_fp'
                        else:  # maccs
                            query_fp = calculate_maccs_fingerprint(query_smiles)
                            fp_col = 'maccs_fp'
                        
                        # Read data
                        db_path = './data.pkl'
                        
                        # Return error if file doesn't exist
                        if not os.path.exists(db_path):
                            return {"error": "Database file does not exist"}

                        db_data = pd.read_pickle(db_path)
                            
                        # Calculate similarity (using Tanimoto coefficient)
                        def calculate_tanimoto(fp):
                            # fp is numpy array
                            # Calculate Tanimoto coefficient
                            intersection = np.sum(np.logical_and(query_fp, fp))
                            union = np.sum(np.logical_or(query_fp, fp))
                            if union == 0:
                                return 0.0
                            return intersection / union
                        
                        db_data['similarity'] = db_data[fp_col].apply(calculate_tanimoto)
                        
                        # Sort by similarity (descending) and Maximum EQE (descending) and return top_n results
                        results = db_data.sort_values('similarity', ascending=False).head(top_n)
                        # Convert maximum_EQE to float type
                        results['maximum_EQE_value'] = results['maximum_EQE_value'].astype(float)
                        results = results.sort_values('maximum_EQE_value', ascending=False).head(1)
                        
                        # Return only the first result
                        return results[['material_SMILES', 'DOI', 'similarity', 'anode', 'hole_injection_layer', 'hole_transport_layer',
                           'emission_layer_details', 'emission_layer_type', 'host', 'dopants',
                           'emission_layers', 'electron_transport_layer',
                           'electron_injection_layer', 'cathode', 'device_emission_wavelength',
                           'device_brightness', 'turn_on_voltage', 'current_efficiency',
                           'power_efficiency', 'maximum_EQE', 'device_lifetime',
                           'device_emission_wavelength_value', 'device_emission_wavelength_unit',
                           'device_brightness_value', 'device_brightness_unit',
                           'turn_on_voltage_value', 'turn_on_voltage_unit',
                           'current_efficiency_value', 'current_efficiency_unit',
                           'power_efficiency_value', 'power_efficiency_unit', 'maximum_EQE_value',
                           'maximum_EQE_unit', 'device_lifetime_value', 'device_lifetime_unit',
                           'pure_emitter', 'dopants_wt_percent', 'dopants_name', 'material_name']]
                
                except ImportError:
                    st.error("RDKit library is required for molecular similarity search. Please run: pip install rdkit")
                    
                # 创建Assistant函数
                def create_assistant(index_id, model_name="qwen-plus"):
                    """Create an Assistant with the specified knowledge base."""
                    assistant = Assistants.create(
                        model=model_name,  # Model list: https://help.aliyun.com/zh/model-studio/getting-started/models
                        name='Organic Light-Emitting Diode (OLED) Laboratory Assistant',
                        description='An assistant for OLED laboratory preparation',
                        instructions='You are an OLED laboratory preparation assistant, specializing in device fabrication and answering all questions about OLED preparation. Use the provided knowledge base to answer user questions. The following information may be helpful: ${documents}. When users provide SMILES strings, you can use the molecular similarity search function to find similar molecules.',
                        tools=[
                            {
                                "type": "rag",  # Specify RAG (Retrieval Augmented Generation) mode
                                "prompt_ra": {
                                    "pipeline_id": [index_id],  # Specify knowledge base index ID
                                    "multiknowledge_rerank_top_n": 10,  # Top N results for multi-knowledge reranking
                                    "rerank_top_n": 5,  # Final top N results after reranking
                                    "parameters": {
                                        "type": "object",
                                        "properties": {
                                            "query_word": {
                                                "type": "str",
                                                "value": "${documents}"  # Dynamic placeholder for query content
                                            }
                                        }
                                    }
                                }
                            },
                            {
                                "type": "function",
                                "function": {
                                    "name": "search_similar_molecules",
                                    "description": "Search for similar molecules based on SMILES string",
                                    "parameters": {
                                        "type": "object",
                                        "properties": {
                                            "smiles": {
                                                "type": "string",
                                                "description": "SMILES string of the molecule"
                                            },
                                            "fp_type": {
                                                "type": "string",
                                                "enum": ["morgan", "maccs"],
                                                "description": "Molecular fingerprint type to use"
                                            },
                                            "top_n": {
                                                "type": "integer",
                                                "description": "Number of similar molecules to return"
                                            }
                                        },
                                        "required": ["smiles"]
                                    }
                                }
                            }
                        ]
                    )
                    return assistant.id
                
                # 初始化会话状态
                if "oled_messages" not in st.session_state:
                    st.session_state.oled_messages = []

                # 初始化 Assistant 和 Thread
                if "assistant_id" not in st.session_state:
                    with st.spinner("Initializing assistant..."):
                        try:
                            assistant_id = create_assistant(st.session_state.knowledge_base_id, model)
                            st.session_state.assistant_id = assistant_id
                        except Exception as e:
                            st.error(f"Failed to initialize assistant: {str(e)}")

                if "thread_id" not in st.session_state:
                    with st.spinner("Creating conversation thread..."):
                        try:
                            thread = Threads.create()
                            st.session_state.thread_id = thread.id
                        except Exception as e:
                            st.error(f"Failed to create conversation thread: {str(e)}")
                
                # 创建聊天界面
                st.subheader("💬 OLED Laboratory Assistant Chat")
                
                # 显示聊天历史
                for message in st.session_state.oled_messages:
                    if message["role"] == "user":
                        with st.chat_message("user"):
                            st.write(message["content"])
                    else:
                        with st.chat_message("assistant"):
                            st.write(message["content"])
                            # 如果有工具调用，显示工具调用信息
                            if "tool_calls" in message and message["tool_calls"]:
                                for tool_call in message["tool_calls"]:
                                    with st.status(f"Tool Call: {tool_call.get('name', 'Knowledge Base Retrieval')}", state="complete"):
                                        st.write("Call Parameters:")
                                        if "args" in tool_call:
                                            st.json(tool_call["args"])
                                        else:
                                            st.write(tool_call.get("content", "No call parameters"))
                                        
                                        if "output" in tool_call:
                                            st.write("Call Results:")
                                            st.write(tool_call["output"])
                
                # 获取用户输入
                if prompt := st.chat_input("Enter your question..."):
                    # 添加用户消息到历史记录
                    user_message = {"role": "user", "content": prompt}
                    st.session_state.oled_messages.append(user_message)
                    with st.chat_message("user"):
                        st.write(prompt)
                    
                    # 显示AI思考中状态
                    with st.chat_message("assistant"):
                        message_placeholder = st.empty()
                        message_placeholder.markdown("Thinking...")
                        
                        # 创建一个容器用于显示检索状态
                        retrieval_status_container = st.container()
                        
                        try:
                            # 创建用户消息
                            Messages.create(thread_id=st.session_state.thread_id, content=prompt)
                            
                            # 使用原生流式输出创建运行
                            run = Runs.create(
                                thread_id=st.session_state.thread_id, 
                                assistant_id=st.session_state.assistant_id,
                                stream=True  # Enable streaming output
                            )

                            # 初始化变量
                            full_response = ""
                            tool_calls_info = []
                            is_tool_call = False
                            tool_call_statuses = {}  # Store tool call status components

                            # 使用 while 循环来处理可能的多轮工具调用
                            while True:
                                # 处理流式输出
                                for event, data in run:
                                    # 处理消息增量更新
                                    if event == 'thread.message.delta':
                                        if hasattr(data, 'delta') and hasattr(data.delta, 'content'):
                                            if hasattr(data.delta.content, 'text') and hasattr(data.delta.content.text, 'value'):
                                                # Update message content
                                                text_delta = data.delta.content.text.value
                                                full_response += text_delta
                                                message_placeholder.markdown(full_response + "▌")
                                    
                                    # 处理工具调用
                                    elif event == 'thread.run.requires_action':
                                        is_tool_call = True
                                        
                                        # Collect tool call information
                                        tool_outputs = []
                                        for tool_call in data.required_action.submit_tool_outputs.tool_calls:
                                            tool_call_id = tool_call.id
                                            
                                            # Create tool call status
                                            with retrieval_status_container:
                                                tool_call_name = tool_call.function.name if hasattr(tool_call, 'function') else "Knowledge Base Retrieval"
                                                tool_call_statuses[tool_call_id] = st.status(
                                                    f"Tool Call: {tool_call_name}", 
                                                    state="running"
                                                )
                                                
                                                # Display call parameters
                                                with tool_call_statuses[tool_call_id]:
                                                    st.write("Call Parameters:")
                                                    try:
                                                        args_json = json.loads(tool_call.function.arguments)
                                                        st.json(args_json)
                                                    except:
                                                        st.write(tool_call.function.arguments)
                                            
                                            # Collect tool call information
                                            tool_call_info = {
                                                "id": tool_call_id,
                                                "type": tool_call.type,
                                                "name": tool_call.function.name if hasattr(tool_call, 'function') else "Knowledge Base Retrieval",
                                                "args": tool_call.function.arguments if hasattr(tool_call, 'function') else "",
                                                "content": tool_call.function.arguments if hasattr(tool_call, 'function') else ""
                                            }
                                            tool_calls_info.append(tool_call_info)
                                            
                                            # Process different types of tool calls
                                            if hasattr(tool_call, 'function') and tool_call.function.name == "search_similar_molecules":
                                                # Process molecular similarity search
                                                args = json.loads(tool_call.function.arguments)
                                                smiles = args.get("smiles", "")
                                                fp_type_arg = args.get("fp_type", fp_type)  # Use parameter or default value
                                                top_n_arg = args.get("top_n", top_n)  # Use parameter or default value
                                                
                                                results = search_similar_molecules(
                                                    query_smiles=smiles,
                                                    fp_type=fp_type_arg,
                                                    top_n=top_n_arg
                                                )
                                                
                                                # Convert results to appropriate format
                                                if isinstance(results, dict) and "error" in results:
                                                    output = results
                                                elif isinstance(results, dict):
                                                    results_df = pd.DataFrame(results)
                                                    output = {"results": results_df.to_dict(orient='records')}
                                                else:
                                                    # If already a DataFrame, use directly
                                                    output = {"results": results.to_dict(orient='records')}
                                                
                                                tool_call_info["output"] = json.dumps(output)
                                                
                                                tool_outputs.append({
                                                    "tool_call_id": tool_call_id,
                                                    "output": json.dumps(output)
                                                })
                                            else:
                                                # Default to knowledge base retrieval
                                                tool_outputs.append({
                                                    "tool_call_id": tool_call_id,
                                                    "output": "Retrieval successful"
                                                })
                                        
                                        # Update tool call status
                                        for tool_call in tool_calls_info:
                                            if tool_call["id"] in tool_call_statuses:
                                                with tool_call_statuses[tool_call["id"]]:
                                                    st.write("Call Results:")
                                                    if "output" in tool_call:
                                                        try:
                                                            output_json = json.loads(tool_call["output"])
                                                            if "results" in output_json:
                                                                # Display similar molecule results
                                                                results_df = pd.DataFrame(output_json["results"])
                                                                st.dataframe(results_df)
                                                            elif "error" in output_json:
                                                                st.error(output_json["error"])
                                                            else:
                                                                st.write(tool_call["output"])
                                                        except:
                                                            st.write(tool_call["output"])
                                                    else:
                                                        st.write("Call completed")
                                                tool_call_statuses[tool_call["id"]].update(state="complete")
                                        
                                        # Submit tool outputs and get new run object
                                        try:
                                            run = Runs.submit_tool_outputs(
                                                thread_id=st.session_state.thread_id,
                                                run_id=data.id,
                                                tool_outputs=tool_outputs,
                                                stream=True  # Enable streaming output
                                            )
                                            # Break out of current for loop, continue processing event stream with new run object
                                            break
                                        except Exception as e:
                                            st.error(f"Failed to submit tool outputs: {str(e)}")
                                            break
                                    
                                    # 处理运行完成事件
                                    elif event == 'thread.run.completed':
                                        # If no complete response, get final message
                                        if not full_response:
                                            # Get message list
                                            msgs = Messages.list(st.session_state.thread_id)
                                            if msgs and 'data' in msgs and len(msgs['data']) > 0:
                                                final_reply = msgs['data'][0]['content'][0]['text']['value']
                                                message_placeholder.markdown(final_reply)
                                                full_response = final_reply
                                        # Run completed, exit loop
                                        break
                                    
                                    # 处理运行失败事件
                                    elif event in ['thread.run.failed', 'thread.run.cancelled', 'thread.run.expired']:
                                        message_placeholder.error(f"Run failed: {event}")
                                        break
                                else:
                                    # For loop completed normally (no break triggered), all events processed
                                    break
                                
                                # Check if need to continue processing new run object
                                if event != 'thread.run.requires_action':
                                    # If last event is not a tool call, no need to continue processing
                                    break
                            
                            # Remove cursor
                            message_placeholder.markdown(full_response)
                            
                            # Add assistant reply to history
                            assistant_message = {
                                "role": "assistant", 
                                "content": full_response,
                                "tool_calls": tool_calls_info
                            }
                            st.session_state.oled_messages.append(assistant_message)
                            
                        except Exception as e:
                            st.error(f"Error processing request: {str(e)}")
                            import traceback
                            st.error(traceback.format_exc())  # Display detailed error information for debugging
                            st.warning("Please check if your API Key and knowledge base ID are correct, or if the DashScope service is available.")
            
            except ImportError:
                st.error("Please install dashscope library to use OLED assistant features: pip install dashscope")
    
    with tab2:
        st.header("Machine Learning Model Training")
        st.selectbox("Model Type", ["Classification", "Regression", "Clustering"])
        st.number_input("Epochs", 1, value=10)
        st.selectbox("Framework", ["PyTorch", "TensorFlow", "Scikit-learn"])
        st.slider("Learning Rate", 0.0001, 0.1, 0.001)

# ---------------------------------------------------------------------------
def main():
    if "current_page" not in st.session_state:
        st.session_state.current_page = "Home"

    page = sidebar_nav()
    st.session_state.current_page = page

    if page == "Home":
        show_home()
    elif page == "Paper Collection Agent":
        show_paper_download()
    elif page == "Information Extraction Agent":
        show_paper_parser()
    elif page == "Deep Research Agent":
        show_machine_learning()

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
