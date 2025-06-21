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
import zipfile


# 设置密码保护
password = st.text_input("请输入密码（私信转账给小闪电获取密码)", type="password")

# 定义正确的密码
correct_password = "flashhhh"  # 替换为你的密码

# 如果输入的密码正确，则加载应用程序
if password == correct_password:
    # 在这里继续你的应用程序代码
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
    if st.sidebar.button("批量分子对接"):
        st.session_state['page'] = "批量分子对接"
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
        st.markdown("### 这个App用到的人工智能算法：")
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
        
        # 使用选项卡分隔单个处理和批量处理
        tab1, tab2 = st.tabs(["单分子处理", "批量处理"])

        with tab1:
            st.info("请上传一个 SDF 文件，或在画布中绘制分子结构/粘贴SMILES")
            st.markdown("**上传**分子文件（SDF 格式）：")
            sdf_file = st.file_uploader("单分子SDF上传", type=["sdf"], key="single_sdf")
            
            st.markdown("**或者** 在下方绘制分子结构/粘贴SMILES：")
            smiles_input = st_ketcher(key="single_ketcher")

            def process_and_show_mol(
                mol: Chem.Mol, 
                uploaded_sdf_name: str = None, 
                user_defined_filename: str = None
            ):
                """处理并展示单个分子"""
                if not mol:
                    return

                try:
                    # 处理手性
                    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                    
                    # 2D 可视化
                    st.subheader("2D 分子结构")
                    img = Draw.MolToImage(mol, size=(300, 300))
                    st.image(img, use_container_width=False)

                    # 生成3D结构
                    mol_3d = Chem.AddHs(mol)
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 42
                    
                    # 检测手性中心
                    chiral_centers = Chem.FindMolChiralCenters(mol_3d)
                    params.enforceChirality = bool(chiral_centers)
                    
                    # 构象生成与优化
                    embed_result = AllChem.EmbedMolecule(mol_3d, params)
                    if embed_result == -1:
                        raise ValueError("构象生成失败")
                    AllChem.MMFFOptimizeMolecule(mol_3d)

                    # 3D 可视化
                    st.subheader("3D 分子结构")
                    mol_block = Chem.MolToMolBlock(mol_3d)
                    xyzview = py3Dmol.view(width=500, height=400)
                    xyzview.addModel(mol_block, "mol")
                    xyzview.setStyle({'stick': {}})
                    xyzview.zoomTo()
                    showmol(xyzview, height=400, width=500)

                    # 生成下载文件
                    if uploaded_sdf_name:
                        base_name = os.path.splitext(uploaded_sdf_name)[0]
                        out_filename = f"{base_name}_prepared.sdf"
                    else:
                        clean_name = re.sub(r'[^a-zA-Z0-9_]', '', user_defined_filename.strip())
                        out_filename = f"{clean_name}.sdf" if user_defined_filename else "ligand_3d.sdf"

                    st.download_button(
                        label="下载3D分子的SDF文件",
                        data=mol_block,
                        file_name=out_filename,
                        mime="chemical/x-mdl-sdfile",
                        use_container_width=True
                    )

                except Exception as e:
                    st.error(f"分子处理失败: {str(e)}")

            # 单个分子处理逻辑
            if sdf_file is not None:
                try:
                    with st.spinner("正在解析SDF文件..."):
                        suppl = Chem.ForwardSDMolSupplier(sdf_file)
                        mol = next((m for m in suppl if m is not None), None)
                        if mol:
                            process_and_show_mol(mol, sdf_file.name)
                        else:
                            st.error("SDF文件中未找到有效分子")
                except Exception as e:
                    st.error(f"SDF文件解析错误: {str(e)}")
            elif smiles_input:
                with st.spinner("正在处理SMILES..."):
                    mol = Chem.MolFromSmiles(smiles_input)
                    if mol:
                        user_name = st.text_input("保存文件名（不含扩展名）", value="my_molecule")
                        process_and_show_mol(mol, user_defined_filename=user_name)
                    else:
                        st.error("无效的SMILES格式")

        with tab2:
            st.subheader("批量处理SMILES文件")
            st.markdown("""
            **CSV文件格式要求：**
            - 必须包含 `mol_name` 和 `smiles` 两列
            - `mol_name` 将用作生成的文件名（仅允许字母、数字和下划线）
            - 示例格式：
            ```
            mol_name | smiles
            ethanol  | CCO
            caffeine | CN1C=NC2=C1C(=O)N(C(=O)N2C)C
            """
            )

            csv_file = st.file_uploader("上传CSV文件", type=["csv"], key="batch_csv")
            
            if csv_file:
                with st.spinner("正在处理批量文件..."):
                    try:
                        df = pd.read_csv(csv_file)
                        if not {'smiles', 'mol_name'}.issubset(df.columns):
                            st.error("CSV文件必须包含'smiles'和'mol_name'两列")
                            st.stop()

                        # 创建临时工作区
                        with tempfile.TemporaryDirectory() as tmpdir:
                            zip_path = os.path.join(tmpdir, "processed_molecules.zip")
                            error_log = []
                            success_count = 0

                            progress_bar = st.progress(0)
                            total = len(df)
                            
                            with zipfile.ZipFile(zip_path, 'w') as zipf:
                                for idx, row in df.iterrows():
                                    try:
                                        # 检查进度
                                        if idx % 1 == 0:
                                            progress = (idx + 1) / total
                                            progress_bar.progress(min(progress, 1.0))
                                        
                                        # 验证分子名称
                                        mol_name = re.sub(r'[^a-zA-Z0-9_]', '', str(row['mol_name']))
                                        if not mol_name:
                                            raise ValueError("无效的分子名称")

                                        # 解析SMILES
                                        mol = Chem.MolFromSmiles(row['smiles'])
                                        if not mol:
                                            raise ValueError("无效的SMILES格式")
                                        
                                        # 处理手性
                                        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                                        chiral_centers = Chem.FindMolChiralCenters(mol)
                                        
                                        # 生成3D结构
                                        mol = Chem.AddHs(mol)
                                        params = AllChem.ETKDGv3()
                                        params.randomSeed = 42
                                        params.enforceChirality = bool(chiral_centers)
                                        
                                        embed_result = AllChem.EmbedMolecule(mol, params)
                                        if embed_result == -1:
                                            raise ValueError("构象生成失败")
                                        AllChem.MMFFOptimizeMolecule(mol)

                                        # 保存文件
                                        sdf_path = os.path.join(tmpdir, f"{mol_name}.sdf")
                                        with Chem.SDWriter(sdf_path) as writer:
                                            writer.write(mol)
                                        
                                        zipf.write(sdf_path, arcname=f"{mol_name}.sdf")
                                        success_count += 1
                                    
                                    except Exception as e:
                                        error_log.append(f"行 {idx+1}: {str(e)}")
                                        continue

                            # 处理完成
                            progress_bar.progress(1.0)
                            
                            # 显示处理结果
                            st.markdown(f"**处理完成** ✅ 成功: {success_count} ❌ 失败: {len(error_log)}")
                            
                            if error_log:
                                with st.expander("错误详情", expanded=False):
                                    st.write("\n".join(error_log))

                            # 提供下载
                            with open(zip_path, "rb") as f:
                                st.download_button(
                                    label="下载处理结果(ZIP)",
                                    data=f,
                                    file_name="processed_molecules.zip",
                                    mime="application/zip",
                                    key="batch_download"
                                )

                    except pd.errors.EmptyDataError:
                        st.error("上传的CSV文件为空")
                    except pd.errors.ParserError:
                        st.error("CSV文件解析失败，请检查文件格式")
                    except Exception as e:
                        st.error(f"文件处理错误: {str(e)}")


    # ------------------------------------------------------------------------------
    # 口袋预测
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # 口袋预测
    # ------------------------------------------------------------------------------


    elif page == "口袋预测":
        import streamlit as st
        import pandas as pd
        import tempfile
        import os
        import hashlib
        from pathlib import Path
        import sh

        TMP_ROOT = Path(tempfile.gettempdir()) / "streamlit-molstar"
        TMP_ROOT.mkdir(exist_ok=True)

        # Helper functions
        def _get_file_type(file_path):
            return os.path.splitext(file_path)[1][1:].lower()

        def get_workspace_info(md5hash, ftype, fresh=False, create=True):
            workdir_p = TMP_ROOT / md5hash
            if fresh and workdir_p.exists():
                sh.rm('-r', str(workdir_p))
            if create:
                workdir_p.mkdir(parents=True, exist_ok=True)
            protein_name = f'{md5hash}-protein.{ftype}'
            protein_file_path_p = workdir_p / protein_name
            pockets_file_path_p = workdir_p / f'{protein_name}_predictions.csv'
            return {
                "workdir": str(workdir_p),
                "protein_file_path": str(protein_file_path_p),
                "pockets_file_path": str(pockets_file_path_p),
            }

        def batch_predict_pockets(protein_file_paths, original_filenames, p2rank_home):
            all_pockets_data = []

            for protein_file_path, original_name in zip(protein_file_paths, original_filenames):
                with open(protein_file_path, 'rb') as f:
                    content = f.read()
                    md5hash = hashlib.md5(content).hexdigest()
                
                ftype = _get_file_type(protein_file_path)
                workspace_info = get_workspace_info(md5hash, ftype, fresh=True, create=True)

                # Copy protein file to workspace
                sh.cp(protein_file_path, workspace_info['protein_file_path'])

                # Run p2rank prediction
                predict_path_p = Path(workspace_info['workdir']) / 'predict'
                predict_path_p.mkdir(parents=True, exist_ok=True)
                cmd = sh.Command(os.path.join(p2rank_home, 'prank'))
                args = ['predict', '-f', workspace_info['protein_file_path'], '-o', str(predict_path_p)]
                cmd(*args, _cwd=p2rank_home, _fg=True)

                protein_file_name = os.path.basename(workspace_info['protein_file_path'])
                tmp_pockets_file_path_p = predict_path_p / f'{protein_file_name}_predictions.csv'
                sh.cp(str(tmp_pockets_file_path_p), workspace_info['pockets_file_path'])

                # Load predictions and use original filename
                df = pd.read_csv(workspace_info['pockets_file_path'])
                df['Protein File'] = original_name  # Use original filename here
                all_pockets_data.append(df)

            if all_pockets_data:
                result_df = pd.concat(all_pockets_data, ignore_index=True)
                return result_df
        st.title("口袋预测")

        # 让用户选择如何加载蛋白质
        option = st.radio("Select how to load the protein:", ( "单个蛋白质口袋预测", "批量蛋白口袋预测", "加载示例蛋白"))

        # 用于保存用户上传的蛋白文件名称（用于替换 Pocket Name）
        uploaded_pdb_filename = None

        if option == "单个蛋白质口袋预测":
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

        elif option == "加载示例蛋白":
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

        elif option == "批量蛋白口袋预测":
            st.header("批量蛋白口袋预测")
            uploaded_files = st.file_uploader("上传多个蛋白质文件（.pdb）", type=["pdb"], accept_multiple_files=True)

            if uploaded_files:
                protein_file_paths = []
                original_filenames = []

                for file in uploaded_files:
                    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                        tmp.write(file.getvalue())
                        tmp.flush()
                        protein_file_paths.append(tmp.name)
                        original_filenames.append(file.name)  # 存储原始文件名

                if st.button("开始批量预测"):
                    progress_bar = st.progress(0)
                    status_text = st.empty()

                    results = []
                    total_files = len(protein_file_paths)

                    for i, (file_path, original_name) in enumerate(zip(protein_file_paths, original_filenames)):
                        status_text.text(f"正在预测：{original_name} ({i+1}/{total_files})")
                        single_result = batch_predict_pockets([file_path], [original_name], p2rank_home='./others/p2rank_2.5/')
                        results.append(single_result)
                        progress_bar.progress((i + 1) / total_files)

                    result_df = pd.concat(results, ignore_index=True)
                    status_text.text("批量预测完成！")
                    st.success("批量预测完成！")

                    # 显示和下载结果
                    st.dataframe(result_df)
                    csv_data = result_df.to_csv(index=False)
                    st.download_button(
                        label="下载所有蛋白预测口袋信息",
                        data=csv_data,
                        file_name="batch_pocket_predictions.csv",
                        mime="text/csv",
                        use_container_width=True
                    )

                # 清理临时文件
                for path in protein_file_paths:
                    if os.path.exists(path):
                        os.remove(path)

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

        size_x = st.number_input("Size X", value=100.0)
        size_y = st.number_input("Size Y", value=100.0)
        size_z = st.number_input("Size Z", value=100.0)

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
    # ------------------------------------------------------------------------------
    # 批量口袋预测与对接
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # 批量口袋预测与对接
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # 批量分子对接（基于已有口袋预测结果）
    # ------------------------------------------------------------------------------
    elif page == "批量分子对接":
        import os
        import uuid
        import time
        import pandas as pd
        import subprocess
        import json
        import zipfile
        from pathlib import Path
        import streamlit as st
        import shutil
        from threading import Thread

        # 初始化任务存储目录
        JOBS_DIR = Path("./jobs")
        JOBS_DIR.mkdir(exist_ok=True)

        # 任务管理函数
        def init_job(job_id):
            """初始化任务目录结构"""
            job_dir = JOBS_DIR / job_id
            (job_dir / "uploads").mkdir(parents=True, exist_ok=True)
            (job_dir / "results").mkdir(exist_ok=True)
            return job_dir

        def save_job_status(job_id, status):
            """保存任务状态"""
            status_file = JOBS_DIR / job_id / "status.json"
            with open(status_file, "w") as f:
                json.dump(status, f)

        def get_job_status(job_id):
            """获取任务状态"""
            status_file = JOBS_DIR / job_id / "status.json"
            try:
                with open(status_file) as f:
                    return json.load(f)
            except:
                return {"state": "not_found", "progress": 0, "log": []}

        # 对接处理线程
        def process_docking_tasks(job_id, tasks_to_run, pockets_df):
            job_dir = JOBS_DIR / job_id
            status = {
                "state": "running",
                "progress": 0,
                "start_time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "log": []
            }
            
            try:
                save_job_status(job_id, status)
                
                # 准备zip文件
                zip_path = job_dir / "results" / "docking_results.zip"
                with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    
                    total_tasks = len(tasks_to_run)
                    for i, task in enumerate(tasks_to_run.iterrows()):
                        idx = i + 1
                        _, row = task
                        protein_path = job_dir / "uploads" / row["Protein"]
                        ligand_path = job_dir / "uploads" / row["Ligand"]
                        
                        # 更新状态
                        status["log"].append(f"任务 {idx}/{total_tasks}: 对接 {row['Protein']} 和 {row['Ligand']}")
                        status["progress"] = idx / total_tasks
                        save_job_status(job_id, status)
                        
                        # 从CSV中获取口袋中心坐标
                        try:
                            # 标准化列名（不区分大小写）
                            protein_col = next(col for col in pockets_df.columns if col.lower() == "protein file")
                            rank_col = next(col for col in pockets_df.columns if col.lower() == "rank")
                            center_x_col = next(col for col in pockets_df.columns if col.lower() == "center_x")
                            center_y_col = next(col for col in pockets_df.columns if col.lower() == "center_y")
                            center_z_col = next(col for col in pockets_df.columns if col.lower() == "center_z")
                            
                            # 获取当前蛋白的rank1口袋
                            protein_filename = row["Protein"]
                            protein_pockets = pockets_df[pockets_df[protein_col] == protein_filename]
                            
                            if len(protein_pockets) == 0:
                                status["log"].append(f"任务 {idx} 错误: 蛋白 {protein_filename} 无口袋数据")
                                continue
                                
                            # 获取rank1口袋
                            rank1_pocket = protein_pockets[protein_pockets[rank_col] == 1].iloc[0]
                            center_coords = [
                                float(rank1_pocket[center_x_col]),
                                float(rank1_pocket[center_y_col]),
                                float(rank1_pocket[center_z_col])
                            ]
                            status["log"].append(f"任务 {idx} 使用口袋中心坐标: {center_coords}")
                            
                        except Exception as e:
                            status["log"].append(f"任务 {idx} 获取口袋坐标失败: {str(e)}")
                            continue

                        # 生成对接网格
                        docking_grid = {
                            "center_x": center_coords[0],
                            "center_y": center_coords[1],
                            "center_z": center_coords[2],
                            "size_x": 25.0,
                            "size_y": 25.0,
                            "size_z": 25.0,
                        }
                        grid_path = job_dir / f"grid_{idx}.json"
                        with open(grid_path, "w") as f:
                            json.dump(docking_grid, f)
                        
                        # 执行对接命令
                        output_ligand_name = f"{Path(ligand_path).stem}_predict"
                        command = (
                            f"python ./others/Uni-Mol/unimol_docking_v2/interface/demo.py "
                            f"--mode single "
                            f"--conf-size 10 "
                            f"--cluster "
                            f"--input-protein {protein_path} "
                            f"--input-ligand {ligand_path} "
                            f"--input-docking-grid {grid_path} "
                            f"--output-ligand-name {output_ligand_name} "
                            f"--output-ligand-dir {job_dir / 'results'} "
                            f"--steric-clash-fix "
                            f"--model-dir ./others/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt"
                        )
                        
                        result = subprocess.run(command, shell=True, capture_output=True, text=True)
                        
                        if result.returncode == 0:
                            # 保存结果到zip
                            output_path = job_dir / "results" / f"{output_ligand_name}.sdf"
                            zip_name = f"{Path(protein_path).stem}_{Path(ligand_path).stem}_docked.sdf"
                            zipf.write(output_path, zip_name)
                            status["log"].append(f"任务 {idx} 对接成功")
                        else:
                            status["log"].append(f"对接失败: {result.stderr}")
                        
                        # 清理临时文件
                        if grid_path.exists():
                            grid_path.unlink()
                
                # 最终状态
                status["state"] = "completed"
                status["end_time"] = time.strftime("%Y-%m-%d %H:%M:%S")
                save_job_status(job_id, status)
                
            except Exception as e:
                status["state"] = "failed"
                status["log"].append(f"处理失败: {str(e)}")
                save_job_status(job_id, status)

        # 页面布局
        st.title("批量分子对接")

        # 任务查询模块
        st.subheader("任务查询")
        job_id = st.text_input("输入任务ID查询状态")
        if job_id:
            status = get_job_status(job_id)
            if status["state"] == "not_found":
                st.error("任务ID不存在")
            else:
                st.write(f"**任务状态**: {status['state'].upper()}")
                st.download_button(
                    "下载结果包",
                    data=(JOBS_DIR / job_id / "results" / "docking_results.zip").read_bytes(),
                    file_name=f"results_{job_id}.zip",
                    disabled=status["state"] != "completed"
                )
                st.subheader("任务日志")
                st.write("\n".join(status["log"]))

        # 新任务提交模块
        st.subheader("新建任务")
        
        # 上传口袋预测结果CSV
        pockets_csv = st.file_uploader(
            "上传批量口袋预测结果CSV文件（请提前在口袋预测模块预测所有蛋白口袋）",
            type=["csv"],
            help="请上传之前批量蛋白口袋预测生成的CSV文件"
        )
        
        if pockets_csv:
            try:
                pockets_df = pd.read_csv(pockets_csv)
                # 标准化列名（不区分大小写）
                pockets_df.columns = [col.strip().lower() for col in pockets_df.columns]
                st.success(f"成功加载口袋数据，共 {len(pockets_df)} 条记录")
                
                # 显示口袋数据预览
                with st.expander("查看口袋数据"):
                    st.dataframe(pockets_df)
                    
                # 上传蛋白质和配体文件
                uploaded_files = st.file_uploader(
                    "上传蛋白质 (PDB) 和配体 (SDF) 文件",
                    type=["pdb", "sdf"],
                    accept_multiple_files=True,
                    key="mol_files"
                )
                
                if uploaded_files:
                    # 生成任务ID
                    job_id = str(uuid.uuid4())[:8]
                    job_dir = init_job(job_id)
                    
                    # 保存上传文件
                    upload_dir = job_dir / "uploads"
                    for f in uploaded_files:
                        (upload_dir / f.name).write_bytes(f.getbuffer())
                    
                    # 生成任务CSV
                    pdb_files = [f.name for f in upload_dir.glob("*.pdb")]
                    sdf_files = [f.name for f in upload_dir.glob("*.sdf")]
                    tasks = []
                    for p in pdb_files:
                        for s in sdf_files:
                            tasks.append({"Protein": p, "Ligand": s, "Run": "Yes"})
                    task_df = pd.DataFrame(tasks)
                    task_csv = task_df.to_csv(index=False)
                    (job_dir / "tasks.csv").write_text(task_csv)
                    
                    # 显示任务信息
                    st.success(f"任务已创建！请保存您的任务ID: **{job_id}**")
                    st.download_button(
                        "下载任务CSV模板",
                        data=task_csv,
                        file_name="tasks.csv",
                        help="编辑后重新上传以修改需要运行的任务"
                    )


                    # 上传修改后的CSV
                    modified_csv = st.file_uploader("上传修改后的任务CSV", type=["csv"])
                    if modified_csv:
                        # 保存修改后的任务文件
                        (job_dir / "modified_tasks.csv").write_bytes(modified_csv.getvalue())
                        st.success("任务文件已更新")
                        
                        # 解析任务
                        try:
                            tasks_df = pd.read_csv(job_dir / "modified_tasks.csv")
                            tasks_to_run = tasks_df[tasks_df["Run"].str.lower() == "yes"]
                            st.write(f"准备运行 {len(tasks_to_run)} 个任务")

                            # 添加任务数量检查提示
                            if len(tasks_to_run) > 20:
                                st.warning("注意：当前任务数超过20个，请修改任务CSV文件减少任务量")

                            if st.button("开始批量对接"):
                                # 添加正式的任务数量检查
                                if len(tasks_to_run) > 20:
                                    st.error(f"错误：任务数量不能超过20个（当前：{len(tasks_to_run)}）。请修改任务文件后重试。")
                                else:
                                    # 启动处理线程
                                    Thread(target=process_docking_tasks, args=(job_id, tasks_to_run, pockets_df)).start()
                                    st.success(f"任务 {job_id} 已开始后台处理，请保存ID用于查询")


                                
                        except Exception as e:
                            st.error(f"CSV解析失败: {str(e)}")
                            
            except Exception as e:
                st.error(f"口袋CSV文件解析失败: {str(e)}")

        # 样式调整
        st.markdown("""
        <style>
        .stDownloadButton > button {
            width: 100%;
            justify-content: center;
        }
        </style>
        """, unsafe_allow_html=True)


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
                    st.write([file.name for file in batch_docking_dir.iterdir()])

                    # 扫描文件夹中的 SDF 和 PDB 文件
                    sdf_files = [file.name for file in batch_docking_dir.glob("*.sdf")]
                    pdb_files = [file.name for file in batch_docking_dir.glob("*.pdb")]

                    if not pdb_files:
                        st.error("在上传的文件中未找到 PDB 文件。请至少上传一个 PDB 文件。")
                    elif not sdf_files:
                        st.error("在上传的文件中未找到 SDF 文件。请至少上传一个 SDF 文件。")
                    else:
                        # 自动生成任务列表（改进后的文件名解析）
                        tasks = []
                        for sdf_file in sdf_files:
                            try:
                                # 使用 split("_", 1) 分割文件名
                                parts = sdf_file.split("_", 1)
                                if len(parts) < 2:
                                    raise ValueError("文件名缺少下划线分隔符")
                                    
                                receptor_name = parts[0]
                                ligand_part = parts[1]
                                
                                # 清理后缀和扩展名
                                ligand_name = ligand_part.replace("_docked", "").split(".", 1)[0]
                                pdb_file = f"{receptor_name}.pdb"
                                
                                tasks.append({
                                    "Protein": pdb_file,
                                    "Ligand": sdf_file,
                                    "Run": "Yes"
                                })
                            except Exception as e:
                                st.warning(f"文件 {sdf_file} 的命名格式不正确（错误：{str(e)}），跳过此文件。")

                        if not tasks:
                            st.error("未生成任何有效的任务，请检查 SDF 文件的命名格式。")
                        else:
                            # 显示自动生成的任务列表
                            st.subheader("自动匹配的任务列表")
                            st.dataframe(pd.DataFrame(tasks))
                            
                            if st.button("开始批量预测", key="start_batch_predicting"):
                                with st.spinner("正在进行批量亲和力预测，请稍候..."):
                                    try:
                                        final_results = []
                                        sdf_files_list = list(batch_docking_dir.glob("*.sdf"))
                                        pdb_files_list = list(batch_docking_dir.glob("*.pdb"))

                                        progress_bar = st.progress(0)
                                        total_files = len(sdf_files_list)
                                        log_messages = []

                                        for i, sdf_file in enumerate(sdf_files_list):
                                            try:
                                                # 改进的文件名解析
                                                stem = sdf_file.stem
                                                parts = stem.split("_", 1)
                                                if len(parts) < 2:
                                                    raise ValueError("文件名格式错误")
                                                    
                                                receptor_name = parts[0]
                                                ligand_part = parts[1].replace("_docked", "")
                                                
                                                pdb_file_name = f"{receptor_name}.pdb"
                                                pdb_file_path = batch_docking_dir / pdb_file_name

                                                if pdb_file_path.exists():
                                                    log_messages.append(f"正在处理：{pdb_file_name} 和 {sdf_file.name}")
                                                    
                                                    with tempfile.TemporaryDirectory() as tmpdir:
                                                        output_csv_path_tmp = Path(tmpdir) / "temp_result.csv"
                                                        
                                                        cmd = [
                                                            "python",
                                                            "./others/PLANET/pred.py",
                                                            "-p", str(pdb_file_path),
                                                            "-l", str(sdf_file),
                                                            "-m", str(sdf_file),
                                                            "-o", str(output_csv_path_tmp)
                                                        ]
                                                        
                                                        result = subprocess.run(cmd, capture_output=True, text=True)
                                                        
                                                        if result.returncode == 0 and output_csv_path_tmp.exists():
                                                            temp_df = pd.read_csv(output_csv_path_tmp)
                                                            if "Binding_Affinity" in temp_df.columns:
                                                                binding_affinity = temp_df["Binding_Affinity"].iloc[0]
                                                                final_results.append({
                                                                    "Protein_File": receptor_name,
                                                                    "Ligand_File": ligand_part,
                                                                    "Binding_Affinity": binding_affinity
                                                                })
                                                                log_messages.append(f"成功：{pdb_file_name} 和 {sdf_file.name}")
                                                            else:
                                                                log_messages.append(f"错误：{sdf_file.name} 结果文件缺少 Binding_Affinity 列")
                                                        else:
                                                            log_messages.append(f"失败：{sdf_file.name}\n错误信息：{result.stderr}")
                                                else:
                                                    log_messages.append(f"跳过：未找到对应的 {pdb_file_name}")

                                                progress_bar.progress((i+1)/total_files)
                                                
                                            except Exception as e:
                                                log_messages.append(f"处理 {sdf_file.name} 时发生错误：{str(e)}")

                                        if final_results:
                                            results_df = pd.DataFrame(final_results)
                                            
                                            # 保存结果
                                            binding_affinity_dir = Path("./Result/Binding_Affinity")
                                            binding_affinity_dir.mkdir(parents=True, exist_ok=True)
                                            
                                            output_csv_path = binding_affinity_dir / "batch_prediction_results.csv"
                                            results_df.to_csv(output_csv_path, index=False)
                                            
                                            # 生成热图
                                            try:
                                                heatmap_data = results_df.pivot(
                                                    index="Protein_File", 
                                                    columns="Ligand_File", 
                                                    values="Binding_Affinity"
                                                )
                                                
                                                plt.figure(figsize=(10, 8), dpi=600)
                                                sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".2f")
                                                plt.title("Binding Affinity Heatmap")
                                                plt.xlabel("Ligands")
                                                plt.ylabel("Proteins")
                                                
                                                heatmap_path = binding_affinity_dir / "binding_affinity_heatmap.png"
                                                plt.savefig(heatmap_path, bbox_inches='tight', dpi=600)
                                                st.pyplot(plt)
                                                
                                                # 更新会话状态
                                                st.session_state.is_batch_predict_done = True
                                                st.session_state.batch_csv_path = str(output_csv_path)
                                                st.session_state.heatmap_path = str(heatmap_path)
                                                st.success("预测完成！")
                                                
                                            except Exception as e:
                                                st.error(f"生成热图失败：{str(e)}")
                                        else:
                                            st.error("未生成有效预测结果")
                                        
                                        # 显示日志
                                        st.subheader("处理日志")
                                        st.text_area("日志详情", value="\n".join(log_messages), height=300)

                                    except Exception as e:
                                        st.error(f"批量预测过程发生严重错误：{str(e)}")

            # 下载结果部分保持不变
            if st.session_state.get('is_batch_predict_done'):
                try:
                    with open(st.session_state.batch_csv_path, "rb") as f_csv, \
                        open(st.session_state.heatmap_path, "rb") as f_img:
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.download_button(
                                label="下载CSV结果",
                                data=f_csv,
                                file_name="batch_predictions.csv",
                                mime="text/csv"
                            )
                        with col2:
                            st.download_button(
                                label="下载热图",
                                data=f_img,
                                file_name="affinity_heatmap.png",
                                mime="image/png"
                            )
                except Exception as e:
                    st.error(f"下载文件失败：{str(e)}")




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