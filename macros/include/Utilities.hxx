#ifndef Macros_Utilities_hxx
#define Macros_Utilities_hxx

#include "Headers.hxx"

/*
 * Get information from the name
 * - input: name
 * - outputs: dirname, filename, listname, objname
 */
void ParseObjectPath(TString name, TString &dirname, TString &filename, TString &listname, TString &objname) {

    dirname.Clear();
    filename.Clear();
    listname.Clear();
    objname.Clear();

    // find where .root is
    TString dot_root = ".root";
    Int_t js = name.Index(dot_root);

    // split name in two:
    // first half: before & including .root, to extract dirname and filename
    TString first_half = name(0, js + dot_root.Length());
    dirname = first_half(0, first_half.Last('/'));
    filename = first_half(first_half.Last('/') + 1, first_half.Length());

    // second half: after .root, to extract listname and objname
    TString second_half = name(js + dot_root.Length(), name.Length() - js);
    listname = second_half(1, second_half.Last('/') - 1);
    objname = second_half(second_half.Last('/') + 1, second_half.Length());

    /*
    std::cout << ".. Utilities :: ParseObjectPath :: dirname  = " << dirname << std::endl;
    std::cout << ".. Utilities :: ParseObjectPath :: filename = " << filename << std::endl;
    std::cout << ".. Utilities :: ParseObjectPath :: listname = " << listname << std::endl;
    std::cout << ".. Utilities :: ParseObjectPath :: objname = " << objname << std::endl;
    */
}

/*
 * Add a single object (from the path) to this list
 */
void AddSingleObject(TList *list, const char *name) {

    TString dirname, filename, listname, objname;
    ParseObjectPath(name, dirname, filename, listname, objname);

    TFile *input_file = new TFile(dirname + "/" + filename, "READ");
    TList *input_list = (TList *)input_file->Get(listname);

    TObject *obj = input_list->FindObject(objname);

    TString ClassName = obj->ClassName();

    // std::cout << ".. Utilities :: AddSingleObject :: Adding " << name << " .." << std::endl;

    if (ClassName == "TH1F") {
        TH1F *input_hist = (TH1F *)input_list->FindObject(objname);
        TH1F *cloned_hist = (TH1F *)input_hist->Clone();
        cloned_hist->SetDirectory(0);
        list->Add(cloned_hist);
    } else if (ClassName == "TH1D") {
        TH1D *input_hist = (TH1D *)input_list->FindObject(objname);
        TH1D *cloned_hist = (TH1D *)input_hist->Clone();
        cloned_hist->SetDirectory(0);
        list->Add(cloned_hist);
    } else if (ClassName == "TH2F") {
        TH2F *input_hist = (TH2F *)input_list->FindObject(objname);
        TH2F *cloned_hist = (TH2F *)input_hist->Clone();
        cloned_hist->SetDirectory(0);
        list->Add(cloned_hist);
    } else if (ClassName == "TTree") {
        TTree *input_tree = (TTree *)input_list->FindObject(objname);
        TTree *cloned_tree = input_tree->CloneTree();
        cloned_tree->SetDirectory(0);
        list->Add(cloned_tree);
    } else {
        std::cerr << "!! Utilities :: AddSingleObject :: Object " << objname << " is not a TH1F, TH1D, TH2F nor TTree !!" << std::endl;
    }

    input_file->Close();
}

/*
 * Add multiple objects (from the path) to this list
 */
void AddObjects(TList *list, const char *name) {

    TString directory, filename, listname, objname;
    ParseObjectPath(name, directory, filename, listname, objname);

    // case A: no wildcard
    if (!filename.MaybeWildcard() && !directory.MaybeWildcard()) {
        AddSingleObject(list, name);
        return;
    }

    // quickly parse the top directory, too
    TString top_directory = directory(0, directory.Last('/'));
    // printf("AddObjects :: top_directory = %s\n", top_directory.Data());

    // prepare the extract all the directories
    TList dir_list;
    const char *one_dir;

    if (directory.MaybeWildcard()) {
        // printf("AddObjects :: wildcard found on directory path\n");

        void *top_dir = gSystem->OpenDirectory(top_directory);

        if (top_dir) {
            TRegexp re(directory, kTRUE);
            while ((one_dir = gSystem->GetDirEntry(top_dir))) {
                if (!strcmp(one_dir, ".") || !strcmp(one_dir, "..")) {
                    continue;
                }
                dir_list.Add(new TObjString(one_dir));
            }
            gSystem->FreeDirectory(top_dir);

            // sort the directories in alphanumeric order
            dir_list.Sort();
        }
    } else {
        // no wilcard on directory, get single dir
        one_dir = directory(directory.Last('/') + 1, directory.Length()).Data();
        dir_list.Add(new TObjString(one_dir));
    }

    // loop over dirs
    TIter next_dir(&dir_list);
    TObjString *current_dir;
    while ((current_dir = (TObjString *)next_dir())) {

        one_dir = current_dir->GetName();

        if (!filename.MaybeWildcard()) {
            AddSingleObject(list, TString::Format("%s/%s/%s/%s/%s", top_directory.Data(), one_dir, filename.Data(), listname.Data(), objname.Data()));
            continue;
        }

        TList file_list;
        const char *one_file;

        void *dir = gSystem->OpenDirectory(top_directory + "/" + (TString)one_dir);
        // const char *epath = gSystem->ExpandPathName(); // in case you want to add fullpath later

        if (dir) {
            // create a TList to store the filenames (not yet sorted)
            TRegexp re(filename, kTRUE);
            while ((one_file = gSystem->GetDirEntry(dir))) {
                if (!strcmp(one_file, ".") || !strcmp(one_file, "..")) {
                    continue;
                }
                TString s = one_file;
                if (filename != one_file && s.Index(re) == kNPOS) {
                    continue;
                }
                // printf("AddObjects :: >> one_file = %s\n", one_file);
                file_list.Add(new TObjString(one_file));
            }
            gSystem->FreeDirectory(dir);

            // sort the files in alphanumeric order
            file_list.Sort();
        } else {
            // printf("somethin happenud\n");
        }

        // loop over files
        TIter next_file(&file_list);
        TObjString *current_file;
        while ((current_file = (TObjString *)next_file())) {

            one_file = current_file->GetName();
            // printf("AddObjects :: Adding %s/%s/%s/%s/%s\n", top_directory.Data(), one_dir, one_file, listname.Data(), objname.Data());
            AddSingleObject(list, TString::Format("%s/%s/%s/%s/%s", top_directory.Data(), one_dir, one_file, listname.Data(), objname.Data()));
        }  // end of loop over files

        file_list.Delete();
    }  // end of loop over dirs

    dir_list.Delete();
}

#endif
