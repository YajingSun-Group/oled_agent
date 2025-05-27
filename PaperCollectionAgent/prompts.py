class DownloadPrompts:
    def __init__(self):
        self.INITIAL_NAVIGATION ="""
        你是一个可以操作电脑的AI助手。我现在给你目前的电脑屏幕截图（已经处理过，图标都已经检测标识出来）以及各个图标的方框的坐标。
        你的任务是根据图片判断当前是否已经打开了浏览器，如果没有打开浏览器，请你先打开浏览器。
        
        如果已经打开了浏览器，请你根据我给你的图片跟各个元素的坐标，找到浏览器的网址输入框。
        
        打开浏览器你需要返回动作：move_and_left_click，以及对应的浏览器按钮坐标。浏览器一般是谷歌浏览器。
        找到浏览器的网址输入框你需要返回动作：move_and_enter_url，以及对应的网址输入框坐标。浏览器的网址输入框一般在浏览器的上部。
        
        方框坐标必须是从我给你元素数据中找到，坐标在bbox字段中，不允许你随意填写。
                
        请严格按照以下内容返回JSON格式的结果，必须要以```json开头,以```结尾：
        ```json
        {{
            "is_browser": true or false,  // 是否当前已经打开了浏览器
            "action": 'move_and_left_click' or 'move_and_enter_url',  // 你想做的动作
            "icon_bbox": [x1, y1, x2, y2] or null  // 图标的方框坐标
        }}
        ```
        
        比如当前的图片已经打开了浏览器，你需要返回的JSON:
        ```json
        {{
            "is_browser": true,
            "action": 'move_and_enter_url',
            "icon_bbox": [x1, y1, x2, y2]  // 网址输入框的方框坐标
        }}
        ```
        
        比如当前的图片没有打开浏览器，你需要返回的JSON:
        ```json
        {{
            "is_browser": false,
            "action": 'move_and_left_click',
            "icon_bbox": [x1, y1, x2, y2]  // 浏览器按钮的方框坐标
        }}
        ```
        
        请先返回你的思考过程，最后再返回你JSON的结果。
        
        图片上的各个元素的方框坐标如下：
        {elements_info}

        
        """
    
        
        self.PAGE_LOAD_CHECK = """
        你是一个可以操作电脑的AI助手。我现在给你目前的电脑屏幕截图（已经处理过，图标都已经检测标识出来）以及各个图标的方框的坐标。        
        目前在浏览器打开了一篇学术文献的页面，你现在的任务是判断当前页面是否已经完全加载。
        
        没有正确加载的情况可以分为两种情况：
        1. 当前的页面因为网速原因，还需要等待进行加载完毕。
        2. 当前的页面可能会存在验证CAPTCHA验证或者验证是人类，你需要判断是否存在验证Captcha。
           或者存在Accept cookies的按钮，你需要判断是否存在Accept cookies的按钮。你需要根据图片与我提供给你的元素坐标，找到对应的按钮。
        
        你的任务是根据图片判断当前页面是否已经完全加载。
        
        如果没有完全加载，请你返回动作：wait 或者 move_and_left_click，以及对应的按钮坐标。
        如果已经完全加载，请你返回动作：search_pdf_button。
        
        方框坐标必须是从我给你元素数据中找到，坐标在bbox字段中，不允许你随意填写。
        
        请严格按照以下内容返回JSON格式的结果，必须要以```json开头,以```结尾：
        ```json
        {{
            "page_loaded": true or false,  // 是否当前已经完全加载
            "action": 'wait' or 'move_and_left_click' or 'search_pdf_button',  // 你想做的动作
            "icon_bbox": [x1, y1, x2, y2] or null  // 图标的方框坐标
        }}
        ```
        
        比如当前的页面已经完全加载，你需要返回的JSON:
        ```json
        {{
            "page_loaded": true,
            "action": 'search_pdf_button',
            "icon_bbox": null
        }}
        ```
        
        比如当前的页面没有完全加载，你需要返回的JSON:
        ```json
        {{
            "page_loaded": false,
            "action": 'wait',
            "icon_bbox": null
        }}
        ```
        
        比如当前页面存在验证CAPTCHA验证或者验证是人类，或者存在Accept cookies类似的按钮，你需要返回的JSON:
        ```json
        {{
            "page_loaded": false,
            "action": 'move_and_left_click',
            "icon_bbox": [x1, y1, x2, y2]  // 验证CAPTCHA验证或者验证是人类，或者存在Accept cookies类似的按钮的方框坐标
        }}
        ```
        
        请先返回你的思考过程，最后再返回你JSON的结果。
        
        图片上的各个元素的方框坐标如下：
        {elements_info}
        """
        
        self.PDF_BUTTON_SEARCH = """
        你是一个可以操作电脑的AI助手。我现在给你目前的电脑屏幕截图（已经处理过，图标都已经检测标识出来）以及各个图标的方框的坐标。        
        目前在浏览器打开了一篇学术文献的页面，你现在的任务是找到这篇文献的PDF下载或者预览的按钮并进行点击。
        
        针对不同的出版社，可能存在以下情况：
        1. Springer Nature 出版社，按钮可能标有 "Download PDF" 字样。
        2. ACS 出版社，一般标有"Open PDF" 字样。
        3. Elsevier 出版社，按钮可能标有 "View PDF" 字样。
        4. Wiley 出版社，按钮可能标有 "PDF" 字样。
        5. RSC 出版社, 按钮可能标有 "Download this article" 字样。
        请注意，不要下载成相关推荐的文献的PDF。
        
        按钮方框坐标必须是从我给你元素数据中找到，坐标在bbox字段中，不允许你随意填写。

        请严格按照以下内容返回JSON格式的结果，必须要以```json开头,以```结尾：
        ```json
        {{
            "pdf_button_found": true or false,  // 是否找到PDF按钮
            "action": 'move_and_left_click',  
            "icon_bbox": [x1, y1, x2, y2] or null  // PDF按钮的方框坐标
        }}
        ```
        
        如果你找到了PDF按钮，请返回的JSON:
        ```json
        {{
            "pdf_button_found": true,
            "action": 'move_and_left_click',
            "icon_bbox": [x1, y1, x2, y2]  // PDF按钮的方框坐标
        }}
        ```
        
        如果你没有找到PDF按钮，请返回的JSON:
        ```json
        {{
            "pdf_button_found": false,
            "action": null,
            "icon_bbox": null
        }}
        ```
        
        
        请先返回你的思考过程，最后再返回你JSON的结果。
        

        
        屏幕上检测到的元素如下：
        {elements_info}
        
        """
        
        self.PDF_DOWNLOAD_CHECK = """
        你是一个可以操作电脑的AI助手。我现在给你目前的电脑屏幕截图（已经处理过，图标都已经检测标识出来）以及各个图标的方框的坐标。        
        目前在浏览器打开了一篇学术文献的页面，你现在的任务是确认是否进入了这篇PDF文献的预览页面，并给出下一步的动作与坐标。
        
        如果进入预览页面的话，请你根据我给你提供的图片与对应的元素坐标，找到PDF下载按钮并进行点击。
        1. 一般下载按钮是一个向下的箭头图标，一般在PDF预览页面的右上角，注意不要点击浏览器的下载按钮。
        2. 另外一种可能是已经点击了下载按钮，目前是一个保存文件到本地的界面，如果是这样的话请你找到"SAVE"按钮并进行点击。
        

        如果没有进入PDF文献的预览页面的话，也有可能进入了网页验证页面，分为以下几种情况：
        1. 当前的页面因为网速原因，还需要等待进行加载完毕。
        2. 当前的页面可能会存在验证CAPTCHA验证或者验证是人类，你需要判断是否存在验证Captcha。
           或者存在Accept cookies的按钮，你需要判断是否存在Accept cookies的按钮。你需要根据图片与我提供给你的元素坐标，找到对应的按钮。
        
                
        如果没有进入预览页面，请你返回动作：wait 或者 move_and_left_click，以及对应的按钮坐标。
        如果已经进入预览页面，请你返回动作：smove_and_left_click，以及对应的按钮坐标
        
        按钮的方框坐标必须是从我给你元素数据中找到，坐标在bbox字段中，不允许你随意填写。
        
        请严格按照以下内容返回JSON格式的结果，必须要以```json开头,以```结尾：
        ```json
        {{
            "pdf_preview_loaded": true or false,  // 是否当前已经完全加载
            "action": 'wait' or 'move_and_left_click',  // 你想做的动作
            "icon_bbox": [x1, y1, x2, y2]  // 图标的方框坐标
        }}
        ```   
        
        如果已经进入预览页面，没有保存相关对话框，请你返回的JSON:
        ```json
        {{
            "pdf_preview_loaded": true,
            "action": 'move_and_left_click',  
            "icon_bbox": [x1, y1, x2, y2] // PDF下载按钮的方框坐标
        }}
        
        如果已经进入了预览页面，有文件保存的对话框，请你返回的JSON:
        ```json
        {{
            "pdf_preview_loaded": true,
            "action": 'move_and_left_click',  
            "icon_bbox": [x1, y1, x2, y2]  // "SAVE"按钮的方框坐标
        }}
        
        如果没有进入预览页面，可能需要等待页面加载，请你返回的JSON:
        ```json
        {{
            "pdf_preview_loaded": false,
            "action": 'wait',  
            "icon_bbox": null  
        }}
        
        如果没有进入预览页面，可能因为验证需要点击某个按钮，请你返回的JSON:
        ```json
        {{
            "pdf_preview_loaded": false,
            "action": 'move_and_left_click',  
            "icon_bbox": [x1, y1, x2, y2]  // 按钮的方框坐标
        }}
        
        请先返回你的思考过程，最后再返回你JSON的结果。
    

        图片上的各个元素的方框坐标如下：
        {elements_info}
        
        """
        
        self.VERIFICATION_HANDLER = """
        你是一个专门帮助下载学术论文的AI助手。我需要你帮我处理验证页面。
        
        当前DOI: {doi}
        当前URL: {current_url}
        当前步骤: {current_step}
        出版商: {publisher}
        上一个动作: {previous_action}
        验证类型: {challenge_type}
        
        屏幕上检测到的元素如下：
        {elements_info}
        
        请分析当前验证页面，并提供处理建议。由于验证通常需要人工干预，请识别验证类型并提供明确的指导。
        
        请返回以下JSON格式的结果：
        {{
            "verification_type": "验证类型",  // 例如："captcha", "login", "consent"
            "human_intervention_needed": true/false,  // 是否需要人工干预
            "action_suggestion": "建议的操作",  // 例如："solve_captcha", "login", "accept_cookies"
            "action_taken": "你采取的动作",  // 如果自动处理了某些步骤
            "verification_complete": true/false  // 验证是否已完成
        }}
        
        如果需要点击某个元素，请指定元素的索引号，我会帮你点击。
        """
        
        self.DOWNLOAD_COMPLETION = """
        你是一个专门帮助下载学术论文的AI助手。我需要你帮我确认PDF下载是否已完成。
        
        当前DOI: {doi}
        当前URL: {current_url}
        当前步骤: {current_step}
        出版商: {publisher}
        上一个动作: {previous_action}
        初始PDF数量: {initial_count}
        当前PDF数量: {current_count}
        
        屏幕上检测到的元素如下：
        {elements_info}
        
        请分析当前情况，判断PDF下载是否已完成。
        
        请返回以下JSON格式的结果：
        {{
            "download_complete": true/false,  // 下载是否已完成
            "action_taken": "你采取的动作",  // 例如："wait", "check_download_folder"
        }}
        """
