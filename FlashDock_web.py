import streamlit as st
import pandas as pd
from streamlit_molstar import st_molstar, st_molstar_rcsb, st_molstar_remote
# 如需使用口袋预测相关函数
from streamlit_molstar.pocket import (
    select_pocket_from_local_protein,
    # 如果你的项目需要也可以 import select_pocket_from_upload_protein
)
# docking 模块
from streamlit_molstar.docking import st_molstar_docking

import os
import json
import subprocess
import tempfile  # 用于创建临时文件
import re
import tqdm

import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
from stmol import showmol
from streamlit_ketcher import st_ketcher

# 如果没有在 session_state 中记录 page，就初始化一个默认值
if 'page' not in st.session_state:
    st.session_state['page'] = '主页'

# 在侧边栏使用按钮来切换页面
st.sidebar.title("Navigation")
if st.sidebar.button("主页"):
    st.session_state['page'] = "主页"

# 新增“准备配体”按钮（插在“主页”和“口袋预测”之间）
if st.sidebar.button("准备配体"):
    st.session_state['page'] = "准备配体"

if st.sidebar.button("口袋预测"):
    st.session_state['page'] = "口袋预测"
if st.sidebar.button("分子对接"):
    st.session_state['page'] = "分子对接"
if st.sidebar.button("批量口袋预测与对接"):
    st.session_state['page'] = "批量口袋预测与对接"
    # 新增“预测亲和力”按钮
if st.sidebar.button("预测亲和力"):
    st.session_state['page'] = "预测亲和力"
# 新增“作者信息”按钮
if st.sidebar.button("作者信息"):
    st.session_state['page'] = "作者信息"

# 获取当前页面
page = st.session_state['page']

# ------------------------------------------------------------------------------
# 主页
# ------------------------------------------------------------------------------
if page == "主页":
    # 使用 HTML 和 Markdown 居中标题
    st.markdown(
        "<h1 style='text-align: center;'>⚡️欢迎使用⚡️</h1>",
        unsafe_allow_html=True,
    )
    st.markdown("<br>", unsafe_allow_html=True)

    # 显示图片 flashdock.png
    if os.path.exists("./others/flashdock.png"):
        st.image("./others/flashdock.png", use_container_width=True)
    else:
        st.error("flashdock.png 文件未找到，请确保它与脚本位于同一目录下。")

    # 在图片和其他内容之间插入若干空行
    st.markdown("<br><br><br><br><br>", unsafe_allow_html=True)

    # 显示 logo.png
    if os.path.exists("./others/logo.png"):
        st.image("./others/logo.png", use_container_width=True)
    else:
        st.error("logo.png 文件未找到，请确保它与脚本位于同一目录下。")

    st.markdown("<br><br><br><br><br><br><br><br>", unsafe_allow_html=True)
    # 在页面底部添加算法信息
    st.markdown("---")  # 分割线
    st.markdown("### 这个App用到的算法：")
    st.markdown("""
    - **对接算法：** [Uni-Mol Docking v2](https://arxiv.org/abs/2405.11769)
    - **口袋预测算法：** [P2Rank](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8)
    - **亲和力预测算法：** [PLANET](https://pubs.acs.org/doi/10.1021/acs.jcim.3c00253)
    """)
    st.markdown("<br><br><br><br>", unsafe_allow_html=True)
    # 下载按钮
    if os.path.exists("./examples/examples.zip"):
        with open("./examples/examples.zip", "rb") as file:
            st.download_button(
                label="下载示例文件",
                data=file,
                file_name="examples.zip",
                mime="application/zip"
            )
    else:
        st.error("examples.zip 文件未找到，请确保它位于 ./examples/ 目录下。")
# ------------------------------------------------------------------------------
# 准备配体
# ------------------------------------------------------------------------------


elif page == "准备配体":
    st.title("准备配体")

    # 1. 允许用户上传一个 SDF 文件
    st.info("请上传一个 SDF 文件，或在画布中绘制分子结构/粘贴SMILES")
    st.markdown("**上传**分子文件（SDF 格式）：")
    sdf_file = st.file_uploader("", type=["sdf"])
    
    # 2. 允许用户使用 Ketcher 绘制或输入 SMILES
    st.markdown("**或者** 在下方绘制分子结构/粘贴SMILES：")
    smiles_input = st_ketcher()

    def process_and_show_mol(
        mol: Chem.Mol, 
        uploaded_sdf_name: str = None, 
        user_defined_filename: str = None
    ):
        """
        对分子进行加氢、3D 嵌入、MMFF 优化并展示 2D/3D 结构；
        根据不同来源决定最终保存的文件名：
        - 如果有 uploaded_sdf_name，则用 "原文件名去除.sdf + '_prepared.sdf'"
        - 如果没有 uploaded_sdf_name，但用户给了自定义文件名，则用 "用户自定义文件名 + '.sdf'"
        """
        if not mol:
            return

        # 2D 可视化
        st.subheader("2D 分子结构")
        st.image(Draw.MolToImage(mol, size=(300, 300)), use_container_width=False)

        # 生成 3D 构象并能量优化
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol_3d)

        # 3D 可视化
        st.subheader("3D 分子结构")
        mol_block = Chem.MolToMolBlock(mol_3d)
        xyzview = py3Dmol.view(width=500, height=400)
        xyzview.addModel(mol_block, "mol")
        xyzview.setStyle({'stick': {}})
        xyzview.zoomTo()
        showmol(xyzview, height=400, width=500)

        # 创建保存为 SDF 格式的临时文件并生成下载链接
        # 创建保存为 SDF 格式的临时文件并生成下载链接
        # 直接提供下载按钮
        if uploaded_sdf_name:
            base_name = os.path.splitext(uploaded_sdf_name)[0]
            out_filename = base_name + "_prepared.sdf"
        else:
            if user_defined_filename:
                out_filename = user_defined_filename.strip() + ".sdf"
            else:
                out_filename = "ligand_3d.sdf"

        # 使用 RDKit 将 3D 分子结构写入 SDF 格式
        sdf_data = Chem.MolToMolBlock(mol_3d)

        # 使用 Streamlit 的 download_button 提供下载功能
        st.download_button(
            label="下载3D分子的SDF文件",
            data=sdf_data,
            file_name=out_filename,
            mime="chemical/x-mdl-sdfile",
            use_container_width=True
        )


    # 首先解析用户上传的 SDF
    mol_from_sdf = None
    uploaded_sdf_name = None

    if sdf_file is not None:
        uploaded_sdf_name = sdf_file.name  # 记录用户上传的文件名
        try:
            sdf_supplier = Chem.ForwardSDMolSupplier(sdf_file)
            mols = [m for m in sdf_supplier if m is not None]
            if len(mols) > 0:
                mol_from_sdf = mols[0]
            else:
                st.error("无法从 SDF 文件中解析出分子，请检查文件格式或内容。")
        except Exception as e:
            st.error(f"读取 SDF 文件出现错误: {e}")

    if mol_from_sdf:
        # 如果成功解析出上传的 SDF，则展示并保存
        process_and_show_mol(mol_from_sdf, uploaded_sdf_name=uploaded_sdf_name)
    else:
        # 如果用户没有上传 SDF 或上传的 SDF 解析失败，则查看 Ketcher 中有没有输入 SMILES
        if smiles_input:
            mol_from_smiles = Chem.MolFromSmiles(smiles_input)
            if mol_from_smiles:
                user_defined_filename = st.text_input("请输入保存时的 SDF 文件名（不含 .sdf）", value="my_mol")
                process_and_show_mol(
                    mol_from_smiles, 
                    uploaded_sdf_name=None, 
                    user_defined_filename=user_defined_filename
                )
            else:
                st.error("SMILES 无效，请重新输入或确认格式。")


# ------------------------------------------------------------------------------
# 口袋预测
# ------------------------------------------------------------------------------
elif page == "口袋预测":
    st.title("口袋预测")

    # 让用户选择如何加载蛋白质
    option = st.radio("Select how to load the protein:", ("上传蛋白质", "加载示例文件"))

    # 用于保存用户上传的蛋白文件名称（用于替换 Pocket Name）
    uploaded_pdb_filename = None

    if option == "上传蛋白质":
        try:
            # 用户上传蛋白质（只出现一次，不会再弹二次上传）
            pdb_file = st.file_uploader("请上传蛋白质文件 (.pdb)", type=["pdb"])
            
            if pdb_file is not None:
                # 记下上传的名称
                uploaded_pdb_filename = pdb_file.name

                # 使用临时文件的方式进行口袋预测
                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                    tmp.write(pdb_file.getvalue())
                    tmp.flush()
                    file_path = tmp.name

                # 调用 p2rank (或其他函数) ，读取该临时文件进行预测
                selected = select_pocket_from_local_protein(
                    file_path,
                    p2rank_home='./others/p2rank_2.5/'
                )
                # 预测完成后删除该临时文件
                os.remove(file_path)

                if selected:
                    pocket = selected
                    st.write("预测到的口袋信息: ", pocket)

                    # 如果 rank=1 的口袋
                    if pocket['rank'] == '1':
                        # 如果上传了文件名，则用之，否则用 pocket['name']
                        final_name = uploaded_pdb_filename if uploaded_pdb_filename else pocket['name']
                        data = {
                            'Pocket Name': [final_name],
                            'Center': [pocket['center']],
                        }
                        df = pd.DataFrame(data)

                        st.write("最优口袋信息预览：")
                        st.dataframe(df)

                        # 用户点击按钮后，才将CSV保存到指定文件夹
                        # 使用 Streamlit 的 download_button 提供单次点击的下载
                        csv_data = df.to_csv(index=False)  # 将 DataFrame 转换为 CSV 数据

                        # 直接显示一个下载按钮
                        st.download_button(
                            label="下载预测口袋信息",
                            data=csv_data,
                            file_name="best_pocket.csv",
                            mime="text/csv",
                            use_container_width=True
                        )

        except Exception as e:
            st.warning(f"处理上传蛋白时发生错误: {e}")

    elif option == "加载示例文件":
        try:
            # 用示例文件名
            uploaded_pdb_filename = "protein_example.pdb"
            # 调用 p2rank 做预测
            selected = select_pocket_from_local_protein(
                "examples/pocket/protein.pdb", 
                p2rank_home='./others/p2rank_2.5/'
            )
            if selected:
                pocket = selected
                st.write("预测到的口袋信息: ", pocket)

                if pocket['rank'] == '1':
                    data = {
                        'Pocket Name': [uploaded_pdb_filename],
                        'Center': [pocket['center']],
                    }
                    df = pd.DataFrame(data)

                    st.write("最优口袋信息预览：")
                    st.dataframe(df)

                    # 将 DataFrame 转换为 CSV 格式的字符串
                    csv_data = df.to_csv(index=False)

                    # 使用 Streamlit 的 download_button 提供单次点击的下载
                    st.download_button(
                        label="下载预测口袋信息",
                        data=csv_data,
                        file_name="best_pocket.csv",
                        mime="text/csv",
                        use_container_width=True
                    )
                        
        except Exception as e:
            st.warning(f"加载示例文件时发生错误: {e}")

# ------------------------------------------------------------------------------
# 分子对接
# ------------------------------------------------------------------------------
elif page == "分子对接":

 

    st.title("分子对接")
    st.write("请上传蛋白质 (PDB 格式) 和配体 (SDF 格式)，并设置对接参数。")

    # 上传蛋白质和配体文件
    protein_file = st.file_uploader("上传蛋白质文件 (.pdb)", type=["pdb"])
    ligand_file = st.file_uploader("上传配体文件 (.sdf)", type=["sdf"])

    # 上传口袋预测结果文件（可选）
    st.write("可选：上传口袋预测结果 CSV 文件，将自动填充对接网格参数。")
    pocket_csv_file = st.file_uploader("上传口袋预测结果文件 (CSV)", type=["csv"])

    # 默认网格参数
    center_x = 0.0
    center_y = 0.0
    center_z = 0.0

    # 处理 CSV 文件获取网格参数
    if pocket_csv_file is not None:
        try:
            pocket_df = pd.read_csv(pocket_csv_file)
            if "Center" in pocket_df.columns:
                center_coords = pocket_df.loc[0, "Center"]
                if isinstance(center_coords, str):
                    coords = [float(c) for c in re.findall(r"[-+]?[0-9]*\.?[0-9]+", center_coords)]
                    if len(coords) == 3:
                        center_x, center_y, center_z = coords
                    else:
                        st.warning("CSV 文件中的 Center 格式不正确，无法自动填充网格参数。")
                else:
                    st.warning("CSV 文件中的 Center 格式不正确，无法自动填充网格参数。")
            else:
                st.warning("CSV 文件中未找到 Center 列，无法自动填充网格参数。")
        except Exception as e:
            st.error(f"读取 CSV 文件时出现错误: {e}")

    # 网格参数输入框
    st.subheader("设置对接口袋参数")
    center_x = st.number_input("Center X", value=center_x)
    center_y = st.number_input("Center Y", value=center_y)
    center_z = st.number_input("Center Z", value=center_z)

    size_x = st.number_input("Size X", value=25.0)
    size_y = st.number_input("Size Y", value=25.0)
    size_z = st.number_input("Size Z", value=25.0)

    # 对接按钮
    if st.button("开始分子对接", key="start_docking"):
        if not protein_file or not ligand_file:
            st.error("请先上传蛋白质 (pdb) 和配体 (sdf) 文件。")
        else:
            try:
                # 创建临时文件夹保存缓存文件
                with tempfile.TemporaryDirectory() as temp_dir:
                    docking_grid = {
                        "center_x": center_x,
                        "center_y": center_y,
                        "center_z": center_z,
                        "size_x": size_x,
                        "size_y": size_y,
                        "size_z": size_z
                    }
                    docking_grid_path = os.path.join(temp_dir, "docking_grid.json")

                    # 保存网格参数为 JSON 文件
                    with open(docking_grid_path, "w") as f:
                        json.dump(docking_grid, f, indent=4)

                    # 保存蛋白质和配体文件到临时目录
                    protein_path = os.path.join(temp_dir, "protein.pdb")
                    ligand_path = os.path.join(temp_dir, "ligand.sdf")

                    with open(protein_path, "wb") as f:
                        f.write(protein_file.getvalue())

                    with open(ligand_path, "wb") as f:
                        f.write(ligand_file.getvalue())

                    # 构造对接命令
                    output_ligand_name = "ligand_predict"
                    command = (
                        f"python ./others/Uni-Mol/unimol_docking_v2/interface/demo.py "
                        f"--mode single "
                        f"--conf-size 10 "
                        f"--cluster "
                        f"--input-protein {protein_path} "
                        f"--input-ligand {ligand_path} "
                        f"--input-docking-grid {docking_grid_path} "
                        f"--output-ligand-name {output_ligand_name} "
                        f"--output-ligand-dir {temp_dir} "
                        f"--steric-clash-fix "
                        f"--model-dir ./others/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt"
                    )

                    # 执行命令
                    result = subprocess.run(command, shell=True, capture_output=True, text=True)

                    # 判断是否成功
                    if result.returncode == 0:
                        st.success("分子对接完成！")
                        st.text_area("对接输出日志", value=result.stdout, height=150)

                        # 处理对接结果文件
                        try:
                            ligand_output_path = os.path.join(temp_dir, f"{output_ligand_name}.sdf")

                            # 重命名输出文件
                            output_name = f"{os.path.splitext(protein_file.name)[0]}_{os.path.splitext(ligand_file.name)[0]}__docked.sdf"
                            renamed_path = os.path.join(temp_dir, output_name)
                            os.rename(ligand_output_path, renamed_path)

                            # 提供下载按钮
                            with open(renamed_path, "rb") as f:
                                sdf_data = f.read()
                                st.download_button(
                                    label="下载对接结果 (SDF)",
                                    data=sdf_data,
                                    file_name=output_name,
                                    mime="chemical/x-mdl-sdfile",
                                    use_container_width=True
                                )

                            # 可视化对接结果
                            st_molstar_docking(
                                protein_path,
                                renamed_path,
                                key="5",
                                height=600
                            )
                        except Exception as e:
                            st.error("处理结果文件时出错，请检查路径或权限。")
                    else:
                        st.error("分子对接失败！")
                        st.text_area("错误信息", value=result.stderr, height=150)

            except Exception as e:
                st.error(f"对接过程出现错误: {e}")


# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
elif page == "批量口袋预测与对接":
    import os
    import pandas as pd
    import subprocess
    import tempfile
    import json
    import zipfile
    from pathlib import Path
    import streamlit as st
    from streamlit_molstar.pocket import select_pocket_from_local_protein
    import shutil

    # 初始化 session_state 变量
    if 'docking_results_zip' not in st.session_state:
        st.session_state['docking_results_zip'] = None
    if 'docked_files' not in st.session_state:
        st.session_state['docked_files'] = []

    st.title("批量口袋预测与分子对接")

    # 用户上传多个文件（pdb或sdf格式）
    uploaded_files = st.file_uploader(
        "上传蛋白质 (PDB 格式) 和配体 (SDF 格式) 文件",
        type=["pdb", "sdf"],
        accept_multiple_files=True
    )

    if uploaded_files:
        # 创建一个缓存文件夹用于保存上传的文件
        with tempfile.TemporaryDirectory() as temp_dir:
            batch_docking_dir = Path(temp_dir)
            os.makedirs(batch_docking_dir, exist_ok=True)

            # 保存上传的文件到临时文件夹
            for uploaded_file in uploaded_files:
                file_path = batch_docking_dir / uploaded_file.name
                with open(file_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())

            # 显示已上传的文件
            st.write(f"已上传文件：")
            st.write([file.name for file in batch_docking_dir.iterdir()])

            # 自动生成任务 CSV 文件
            def generate_task_csv():
                pdb_files = list(batch_docking_dir.glob("*.pdb"))
                sdf_files = list(batch_docking_dir.glob("*.sdf"))

                if not pdb_files:
                    st.error(f"在 {batch_docking_dir} 文件夹中未找到 PDB 文件。请添加至少一个 PDB 文件。")
                    return None
                if not sdf_files:
                    st.error(f"在 {batch_docking_dir} 文件夹中未找到 SDF 文件。请添加至少一个 SDF 文件。")
                    return None

                tasks = []
                for pdb_file in pdb_files:
                    for sdf_file in sdf_files:
                        tasks.append({
                            "Protein": pdb_file.name,
                            "Ligand": sdf_file.name,
                            "Run": "Yes"  # 默认所有任务为 "Yes"
                        })

                task_df = pd.DataFrame(tasks)
                return task_df

            # 生成任务 DataFrame
            task_df = generate_task_csv()

            if task_df is not None:
                # 提供下载任务 CSV 的按钮
                csv = task_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="下载任务 CSV 文件",
                    data=csv,
                    file_name="docking_tasks.csv",
                    mime="text/csv"
                )

                st.markdown("---")
                st.info("""
                    1. 请提前上传蛋白质 (PDB 格式) 和配体 (SDF 格式) 文件。
                    2. 本app会默认上传的所有蛋白和配体进行对接，并生成任务CSV文件。
                    3. 你可以下载并编辑任务CSV文件，修改 `Run` 列为 `No` 的任务将被跳过，`Yes` 的任务将被执行。
                    4. 修改完成后，上传修改后的 CSV 文件并点击“开始批量预测和对接”按钮。
                """)

                # 上传修改后的任务 CSV 文件
                uploaded_csv = st.file_uploader(
                    "上传修改后的任务 CSV 文件",
                    type=["csv"],
                    key="upload_task_csv"
                )

                if uploaded_csv is not None:
                    try:
                        uploaded_tasks_df = pd.read_csv(uploaded_csv)

                        # 检查必要的列是否存在
                        required_columns = {"Protein", "Ligand", "Run"}
                        if not required_columns.issubset(uploaded_tasks_df.columns):
                            st.error(f"上传的任务文件缺少必要的列：{required_columns - set(uploaded_tasks_df.columns)}")
                        else:
                            # 过滤需要运行的任务
                            tasks_to_run = uploaded_tasks_df[uploaded_tasks_df["Run"].str.lower() == "yes"]

                            if tasks_to_run.empty:
                                st.warning("没有任务需要运行，请确保至少有一项任务的 `Run` 列为 `Yes`。")
                            else:
                                st.write(f"发现 {len(tasks_to_run)} 个任务需要运行。")

                                # 显示需要运行的任务表格
                                st.subheader("待运行的任务列表")
                                st.dataframe(tasks_to_run[['Protein', 'Ligand']].reset_index(drop=True))

                                # 开始批量预测和对接按钮
                                if st.button("开始批量预测和对接", key="start_batch_processing"):
                                    log_messages = []
                                    progress_bar = st.progress(0)
                                    status_text = st.empty()

                                    # 临时文件夹用于存储对接结果
                                    with tempfile.TemporaryDirectory() as processing_dir:
                                        # 用于保存对接结果的 zip 文件
                                        zip_file_path = Path(processing_dir) / "docking_results.zip"
                                        with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zipf:

                                            for i, task in tasks_to_run.iterrows():
                                                protein_path = batch_docking_dir / task["Protein"]
                                                ligand_path = batch_docking_dir / task["Ligand"]

                                                # 口袋预测前更新状态
                                                status_text.text(f"任务 {i + 1}/{len(tasks_to_run)}: 正在为 {task['Protein']} 预测口袋...")

                                                # 口袋预测
                                                try:
                                                    pocket_result = select_pocket_from_local_protein(
                                                        str(protein_path),
                                                        p2rank_home='./others/p2rank_2.5/',
                                                        key=f"select_pocket_{i}"
                                                    )
                                                except Exception as e:
                                                    log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 口袋预测失败：{e}")
                                                    progress_bar.progress((i + 1) / len(tasks_to_run))
                                                    continue

                                                if pocket_result:
                                                    try:
                                                        center_coords = [float(coord) for coord in pocket_result["center"]]
                                                    except ValueError:
                                                        log_messages.append(f"任务 {task['Protein']} 的口袋中心坐标格式不正确。")
                                                        progress_bar.progress((i + 1) / len(tasks_to_run))
                                                        continue

                                                    docking_grid = {
                                                        "center_x": center_coords[0],
                                                        "center_y": center_coords[1],
                                                        "center_z": center_coords[2],
                                                        "size_x": 25.0,
                                                        "size_y": 25.0,
                                                        "size_z": 25.0,
                                                    }

                                                    # 创建临时目录存储对接网格
                                                    with tempfile.TemporaryDirectory() as grid_dir:
                                                        docking_grid_path = Path(grid_dir) / "docking_grid.json"
                                                        with open(docking_grid_path, "w") as f:
                                                            json.dump(docking_grid, f, indent=4)

                                                        # 更新状态：开始对接
                                                        status_text.text(f"任务 {i + 1}/{len(tasks_to_run)}: 正在对接 {task['Protein']} 和 {task['Ligand']}...")

                                                        # 构造对接命令
                                                        output_ligand_name = f"{Path(ligand_path).stem}_predict"
                                                        command = (
                                                            f"python ./others/Uni-Mol/unimol_docking_v2/interface/demo.py "
                                                            f"--mode single "
                                                            f"--conf-size 10 "
                                                            f"--cluster "
                                                            f"--input-protein {protein_path} "
                                                            f"--input-ligand {ligand_path} "
                                                            f"--input-docking-grid {docking_grid_path} "
                                                            f"--output-ligand-name {output_ligand_name} "
                                                            f"--output-ligand-dir {processing_dir} "
                                                            f"--steric-clash-fix "
                                                            f"--model-dir ./others/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt"
                                                        )

                                                        # 执行对接命令
                                                        result = subprocess.run(command, shell=True, capture_output=True, text=True)

                                                        if result.returncode == 0:
                                                            ligand_output_path = Path(processing_dir) / f"{output_ligand_name}.sdf"
                                                            output_name = f"{Path(protein_path).stem}_{Path(ligand_path).stem}_docked.sdf"
                                                            renamed_path = Path(processing_dir) / output_name

                                                            try:
                                                                # 将对接结果重命名并添加到 zip 文件中
                                                                shutil.move(str(ligand_output_path), str(renamed_path))
                                                                zipf.write(renamed_path, renamed_path.name)

                                                                # 读取重命名后的文件为字节
                                                                with open(renamed_path, "rb") as f:
                                                                    file_bytes = f.read()

                                                                # 将文件信息添加到 session_state
                                                                st.session_state['docked_files'].append({
                                                                    "name": output_name,
                                                                    "data": file_bytes
                                                                })

                                                                log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 对接完成。结果已添加到 zip 文件。")
                                                            except Exception as e:
                                                                log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 结果保存失败：{e}")
                                                        else:
                                                            log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 对接失败。错误信息：{result.stderr}")

                                                else:
                                                    log_messages.append(f"任务 {task['Protein']} 的口袋信息未找到。")

                                                # 更新进度条
                                                progress_bar.progress((i + 1) / len(tasks_to_run))

                                        # 读取 zip 文件为字节并存储在 session_state
                                        try:
                                            with open(zip_file_path, "rb") as f_zip:
                                                st.session_state['docking_results_zip'] = f_zip.read()
                                        except Exception as e:
                                            st.error(f"读取 ZIP 文件时出错：{e}")

                                    # 显示日志
                                    st.text_area("任务日志", value="\n".join(log_messages), height=300)

                                    st.success("所有任务已完成。你可以下载对接结果。")

                    except pd.errors.EmptyDataError:
                        st.error("上传的 CSV 文件为空，请检查文件内容。")
                    except pd.errors.ParserError:
                        st.error("上传的 CSV 文件格式错误，请确保文件为有效的 CSV 格式。")
                    except Exception as e:
                        st.error(f"读取任务文件时出错：{e}")

    # 显示下载按钮（无论是否有新上传，保持下载按钮可见）
    if st.session_state.get('docking_results_zip') or st.session_state.get('docked_files'):
        st.markdown("---")
        st.subheader("下载对接结果")

        # 下载整个 ZIP 文件
        if st.session_state.get('docking_results_zip'):
            st.download_button(
                label="下载所有对接结果 (ZIP)",
                data=st.session_state['docking_results_zip'],
                file_name="docking_results.zip",
                mime="application/zip"
            )

        # 下载单个对接结果文件
        if st.session_state.get('docked_files'):
            st.markdown("### 单个对接结果下载")
            for docked_file in st.session_state['docked_files']:
                st.download_button(
                    label=f"下载 {docked_file['name']}",
                    data=docked_file['data'],
                    file_name=docked_file['name'],
                    mime="chemical/x-mdl-sdfile",
                    key=f"download_{docked_file['name']}"
                )




# ------------------------------------------------------------------------------
# 预测亲和力
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 预测亲和力
# ------------------------------------------------------------------------------
import os
import tempfile
import subprocess
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# 初始化 session_state
if 'is_batch_predict_done' not in st.session_state:
    st.session_state.is_batch_predict_done = False
    st.session_state.batch_csv_path = ""
    st.session_state.heatmap_path = ""

if page == "预测亲和力":
    st.title("预测亲和力")
    st.write("在此页面，你可以进行小分子与蛋白质的结合亲和力预测。选择单个预测或批量预测模式。")

    # 模式选择
    mode = st.radio("选择模式", ("单个预测", "批量预测"))

    if mode == "单个预测":
        st.subheader("单个蛋白与小分子的亲和力预测")

        # 用户上传蛋白质 PDB 文件
        protein_file = st.file_uploader("上传蛋白质 PDB 文件", type=["pdb"])

        # 用户上传小分子 SDF 文件
        ligand_file = st.file_uploader("上传小分子 SDF 文件", type=["sdf"])

        # 按钮触发预测
        if st.button("开始预测", key="start_predicting"):
            if protein_file is None:
                st.error("请上传蛋白质 PDB 文件。")
            elif ligand_file is None:
                st.error("请上传小分子 SDF 文件。")
            else:
                with st.spinner("正在进行亲和力预测，请稍候..."):
                    try:
                        # 创建临时目录
                        with tempfile.TemporaryDirectory() as tmpdir:
                            # 保存上传的蛋白质文件
                            protein_path = os.path.join(tmpdir, protein_file.name)
                            with open(protein_path, "wb") as f:
                                f.write(protein_file.getbuffer())

                            # 保存上传的小分子文件
                            ligand_path = os.path.join(tmpdir, ligand_file.name)
                            with open(ligand_path, "wb") as f:
                                f.write(ligand_file.getbuffer())

                            # 输出 CSV 文件路径
                            output_csv_path = os.path.join(tmpdir, "single_prediction.csv")

                            # 调用预测脚本
                            pred_dir = "./others/PLANET"
                            pred_script = "pred.py"
                            pred_script_path = os.path.join(pred_dir, pred_script)

                            cmd = [
                                "python",
                                pred_script_path,
                                "-p", protein_path,
                                "-l", ligand_path,
                                "-m", ligand_path,
                                "-o", output_csv_path
                            ]

                            result = subprocess.run(cmd, capture_output=True, text=True)

                            if result.returncode != 0:
                                st.error(f"预测过程中发生错误:\n{result.stderr}")
                            else:
                                if os.path.exists(output_csv_path):
                                    df = pd.read_csv(output_csv_path)
                                    st.success("预测完成！结果如下：")
                                    st.dataframe(df)
                                else:
                                    st.error("预测完成但未找到输出 CSV 文件。")
                    except Exception as e:
                        st.error(f"发生异常: {e}")

    elif mode == "批量预测":
        st.subheader("批量蛋白与小分子亲和力预测")

        # 用户上传多个文件（PDB 和 SDF 格式）
        uploaded_files = st.file_uploader(
            "上传蛋白质 (PDB 格式) 和配体 (SDF 格式) 文件",
            type=["pdb", "sdf"],
            accept_multiple_files=True
        )

        if uploaded_files:
            # 创建临时目录并保存上传的文件
            with tempfile.TemporaryDirectory() as temp_dir:
                batch_docking_dir = Path(temp_dir)
                os.makedirs(batch_docking_dir, exist_ok=True)

                # 保存上传的文件到临时文件夹
                for uploaded_file in uploaded_files:
                    file_path = batch_docking_dir / uploaded_file.name
                    with open(file_path, "wb") as f:
                        f.write(uploaded_file.getbuffer())

                # 显示已上传的文件
                st.write("已上传文件：")
                # 使用 Path.iterdir() 返回 Path 对象，具有 name 属性
                st.write([file.name for file in batch_docking_dir.iterdir()])

                # 扫描文件夹中的 SDF 和 PDB 文件
                sdf_files = [file.name for file in batch_docking_dir.glob("*.sdf")]
                pdb_files = [file.name for file in batch_docking_dir.glob("*.pdb")]

                if not pdb_files:
                    st.error("在上传的文件中未找到 PDB 文件。请至少上传一个 PDB 文件。")
                elif not sdf_files:
                    st.error("在上传的文件中未找到 SDF 文件。请至少上传一个 SDF 文件。")
                else:
                    # 自动生成任务列表
                    tasks = []
                    for sdf_file in sdf_files:
                        try:
                            receptor_name = sdf_file.split("_")[0]
                            ligand_name = sdf_file.split("_")[1].split(".")[0]  # 去除扩展名
                            pdb_file = receptor_name + ".pdb"
                            tasks.append({
                                "Protein": pdb_file,
                                "Ligand": sdf_file,
                                "Run": "Yes"  # 默认所有任务为 "Yes"
                            })
                        except IndexError:
                            st.warning(f"文件 {sdf_file} 的命名格式不正确，跳过此文件。")

                    if not tasks:
                        st.error("未生成任何有效的任务，请检查 SDF 文件的命名格式。")
                    else:
                        # 提供开始批量预测按钮
                        if st.button("开始批量预测", key="start_batch_predicting"):
                            with st.spinner("正在进行批量亲和力预测，请稍候..."):
                                try:
                                    final_results = []

                                    # 扫描文件夹中的 SDF 和 PDB 文件
                                    sdf_files_list = [file for file in batch_docking_dir.glob("*.sdf")]
                                    pdb_files_list = [file for file in batch_docking_dir.glob("*.pdb")]

                                    st.write("发现以下蛋白质文件：")
                                    st.write([file.name for file in pdb_files_list])
                                    st.write("发现以下配体文件：")
                                    st.write([file.name for file in sdf_files_list])

                                    progress_bar = st.progress(0)
                                    total_files = len(sdf_files_list)
                                    current_task = 0  # 用于计算进度

                                    log_messages = []

                                    for i, sdf_file in enumerate(sdf_files_list):
                                        receptor_name = sdf_file.stem.split("_")[0]
                                        ligand_name = sdf_file.stem.split("_")[1]
                                        pdb_file_name = receptor_name + ".pdb"
                                        pdb_file_path = batch_docking_dir / pdb_file_name
                                        sdf_file_path = sdf_file

                                        if pdb_file_path.exists():
                                            st.text(f"正在计算第 {i + 1}/{total_files} 对：蛋白 {pdb_file_path} 和 配体 {sdf_file_path} 的亲和力...")
                                            try:
                                                with tempfile.TemporaryDirectory() as tmpdir:
                                                    output_csv_path_tmp = os.path.join(tmpdir, "temp_result.csv")

                                                    cmd = [
                                                        "python",
                                                        "./others/PLANET/pred.py",
                                                        "-p", str(pdb_file_path),
                                                        "-l", str(sdf_file_path),
                                                        "-m", str(sdf_file_path),
                                                        "-o", output_csv_path_tmp
                                                    ]

                                                    result = subprocess.run(cmd, capture_output=True, text=True)

                                                    if result.returncode == 0 and os.path.exists(output_csv_path_tmp):
                                                        temp_df = pd.read_csv(output_csv_path_tmp)
                                                        if "Binding_Affinity" in temp_df.columns:
                                                            binding_affinity = temp_df["Binding_Affinity"].iloc[0]
                                                            final_results.append({
                                                                "Protein_File": receptor_name,
                                                                "Ligand_File": ligand_name,
                                                                "Binding_Affinity": binding_affinity
                                                            })
                                                            log_messages.append(f"任务 {pdb_file_name} 和 {sdf_file.name} 预测成功。")
                                                        else:
                                                            log_messages.append(f"任务 {pdb_file_name} 和 {sdf_file.name} 预测结果缺少 'Binding_Affinity' 列。")
                                                    else:
                                                        log_messages.append(f"任务 {pdb_file_name} 和 {sdf_file.name} 预测失败。")
                                            except Exception as e:
                                                log_messages.append(f"任务 {pdb_file_name} 和 {sdf_file.name} 预测过程中发生异常：{e}")
                                        else:
                                            log_messages.append(f"任务 {pdb_file_name} 不存在，跳过。")

                                        # 更新进度条
                                        current_task += 1
                                        progress_bar.progress(current_task / total_files)

                                    if final_results:
                                        results_df = pd.DataFrame(final_results)

                                        # 保存结果到 Binding_Affinity 文件夹
                                        binding_affinity_dir = Path("./Result/Binding_Affinity")
                                        binding_affinity_dir.mkdir(parents=True, exist_ok=True)

                                        # 保存批量预测结果 CSV 文件
                                        output_csv_path = binding_affinity_dir / "batch_prediction_results.csv"
                                        results_df.to_csv(output_csv_path, index=False)

                                        # 绘制热图并保存
                                        st.subheader("亲和力热图")
                                        try:
                                            heatmap_data = results_df.pivot(index="Protein_File", columns="Ligand_File", values="Binding_Affinity")
                                            plt.figure(figsize=(10, 8), dpi=600)
                                            sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".2f")
                                            plt.xlabel("Ligands")
                                            plt.ylabel("Proteins")
                                            plt.title("Binding Affinity Heatmap")

                                            # 保存热图到 PNG 文件
                                            heatmap_path = binding_affinity_dir / "binding_affinity_heatmap.png"
                                            plt.savefig(heatmap_path, dpi=600)
                                            st.pyplot(plt)  # This will display the heatmap

                                            # plt.close()

                                            # 在 session_state 中保存文件路径
                                            st.session_state.is_batch_predict_done = True
                                            st.session_state.batch_csv_path = str(output_csv_path)
                                            st.session_state.heatmap_path = str(heatmap_path)

                                            st.success("批量预测完成！你可以下载预测结果和热图。")
                                        except Exception as e:
                                            st.error(f"绘制热图时发生错误：{e}")
                                    else:
                                        st.error("未生成任何预测结果。")

                                    # 显示日志
                                    st.text_area("任务日志", value="\n".join(log_messages), height=300)

                                except Exception as e:
                                    st.error(f"批量预测过程中发生错误：{e}")

        # 如果批量预测完成，显示下载按钮
        if st.session_state.is_batch_predict_done:
            try:
                with open(st.session_state.batch_csv_path, "rb") as f_csv, open(st.session_state.heatmap_path, "rb") as f_img:
                    # 提供批量预测结果和热图的下载按钮
                    st.download_button(
                        label="下载批量预测结果 CSV",
                        data=f_csv,
                        file_name="batch_prediction_results.csv",
                        mime="text/csv"
                    )

                    st.download_button(
                        label="下载亲和力热图 PNG",
                        data=f_img,
                        file_name="binding_affinity_heatmap.png",
                        mime="image/png"
                    )
            except Exception as e:
                st.error(f"读取预测结果文件时发生错误：{e}")




# ------------------------------------------------------------------------------
# 作者信息
# ------------------------------------------------------------------------------
if page == "作者信息":
    # 使用 HTML 和 Markdown 居中标题
    st.markdown(
        "<h1 style='text-align: center;'>👽小闪电-FLASH⚡️</h1>",
        unsafe_allow_html=True,
    )
    st.markdown("<br>", unsafe_allow_html=True)

    # 显示图片 flashdock.png
    if os.path.exists("./others/author.png"):
        st.image("./others/author.png", use_container_width=True)
    else:
        st.error("author.png 文件未找到，请确保它与脚本位于同一目录下。")