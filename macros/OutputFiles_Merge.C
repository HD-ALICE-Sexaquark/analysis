#include "include/Headers.hxx"

void OutputFiles_Merge(TString InputDir, TString OutputFile = "AnalysisResults_merged.root") {

    // get all files that are named AnalysisResults_* within the InputDir
    std::vector<TString> InputFiles;
    TSystemDirectory dir(InputDir, InputDir);
    TList *files = dir.GetListOfFiles();
    TIter next(files);
    TSystemFile *file;
    TString filename;
    while ((file = (TSystemFile *)next())) {
        filename = file->GetName();
        if (filename.Contains("AnalysisResults_")) {
            std::cout << "Adding file: " << filename << std::endl;
            InputFiles.push_back(filename);
        }
    }

    std::cout << "Merging " << InputFiles.size() << " files!" << std::endl;
}
