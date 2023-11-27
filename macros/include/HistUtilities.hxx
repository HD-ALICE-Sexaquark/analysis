#ifndef HistUtilities_hxx
#define HistUtilities_hxx

#include "Headers.hxx"
/**
  Get information from the name
  - input: name
  - outputs: dirname, filename, listname, histname
**/
void ParseHistFilename(TString name, TString &dirname, TString &filename, TString &listname, TString &histname) {

    dirname.Clear();
    filename.Clear();
    listname.Clear();
    histname.Clear();

    // find where .root is
    TString dot_root = ".root";
    Int_t js = name.Index(dot_root);

    // split name in two:
    // first half: before & including .root, to extract dirname and filename
    TString first_half = name(0, js + dot_root.Length());
    dirname = first_half(0, first_half.Last('/'));
    filename = first_half(first_half.Last('/') + 1, first_half.Length());

    // second half: after .root, to extract listname and histname
    TString second_half = name(js + dot_root.Length(), name.Length() - js);
    listname = second_half(1, second_half.Last('/') - 1);
    histname = second_half(second_half.Last('/') + 1, second_half.Length());

    // (debug)
    std::cout << "ParseHistFilename :: dirname  = " << dirname << std::endl;
    std::cout << "ParseHistFilename :: filename = " << filename << std::endl;
    std::cout << "ParseHistFilename :: listname = " << listname << std::endl;
    std::cout << "ParseHistFilename :: histname = " << histname << std::endl;
}

/**
  Add a new file to this list
**/
void AddFileToList(TList *list, const char *name) {
    TString dirname, filename, listname, histname;
    ParseHistFilename(name, dirname, filename, listname, histname);

    TFile *input_file = new TFile(dirname + "/" + filename, "READ");
    TList *input_list = (TList *)input_file->Get(listname);
    TH1 *input_hist = (TH1 *)input_list->FindObject(histname);

    list->Add(input_hist);
}

/**

**/
void AddHistsToList(TList *list, const char *name) {

    TString directory, filename, listname, histname;
    ParseHistFilename(name, directory, filename, listname, histname);

    // case A: no wildcard
    if (!filename.MaybeWildcard() && !directory.MaybeWildcard()) {
        printf("AddHistsToList :: wildcard not found\n");
        printf("AddHistsToList :: Adding %s\n", name);
        AddFileToList(list, name);
        return;
    }

    // quickly parse the top directory, too
    TString top_directory = directory(0, directory.Last('/'));
    // printf("AddTreesToList :: top_directory = %s\n", top_directory.Data());

    // prepare the extract all the directories
    TList dir_list;
    const char *one_dir;

    if (directory.MaybeWildcard()) {
        // printf("AddTreesToList :: wildcard found on directory path\n");

        void *top_dir = gSystem->OpenDirectory(top_directory);

        if (top_dir) {
            TRegexp re(directory, kTRUE);
            // printf("AddTreesToList :: opened top_directory\n");
            while ((one_dir = gSystem->GetDirEntry(top_dir))) {
                // printf("AddTreesToList :: one_dir = %s\n", one_dir);
                if (!strcmp(one_dir, ".") || !strcmp(one_dir, "..")) {
                    // printf("AddTreesToList :: -> excluded 1\n");
                    continue;
                }
                /*
                        TString s = one_dir;
                        if (directory != one_dir && s.Index(re) == kNPOS) {
                          printf("AddTreesToList :: -> excluded 2\n");
                          continue;
                        }
                 */
                // printf("AddTreesToList :: >> one_dir = %s\n", one_dir);
                dir_list.Add(new TObjString(one_dir));
            }
            gSystem->FreeDirectory(top_dir);

            // sort the directories in alphanumeric order
            dir_list.Sort();
        }
    } else {
        // no wilcard on directory, get only di
        one_dir = directory(directory.Last('/') + 1, directory.Length()).Data();
        dir_list.Add(new TObjString(one_dir));
    }

    // empty line
    // printf("AddTreesToList ::\n");

    // loop over dirs
    TIter next_dir(&dir_list);
    TObjString *current_dir;
    while ((current_dir = (TObjString *)next_dir())) {

        one_dir = current_dir->GetName();
        // printf("AddTreesToList :: td+od = %s/%s\n", top_directory.Data(), one_dir);

        // no wildcard on file
        if (!filename.MaybeWildcard()) {
            // printf("AddTreesToList :: wildcard not found on filename\n");
            // printf("AddTreesToList :: Adding %s/%s/%s/%s/%s\n", top_directory.Data(), one_dir, filename.Data(), listname.Data(),
            // histname.Data());
            AddFileToList(
                list, TString::Format("%s/%s/%s/%s/%s", top_directory.Data(), one_dir, filename.Data(), listname.Data(), histname.Data()));
            continue;
        }

        // if there's a wildcard on the file
        // printf("AddTreesToList :: wildcard found on filename\n");
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
                // printf("AddTreesToList :: >> one_file = %s\n", one_file);
                file_list.Add(new TObjString(one_file));
            }
            gSystem->FreeDirectory(dir);

            // sort the files in alphanumeric order
            file_list.Sort();
        } else {
            // printf("somethin happenud\n");
        }

        // empty line
        // printf("AddTreesToList ::\n");

        // loop over files
        TIter next_file(&file_list);
        TObjString *current_file;
        while ((current_file = (TObjString *)next_file())) {

            one_file = current_file->GetName();
            // printf("AddTreesToList :: Adding %s/%s/%s/%s/%s\n", top_directory.Data(), one_dir, one_file, listname.Data(),
            // histname.Data());
            AddFileToList(list,
                          TString::Format("%s/%s/%s/%s/%s", top_directory.Data(), one_dir, one_file, listname.Data(), histname.Data()));
        }  // end of loop over files

        file_list.Delete();
    }  // end of loop over dirs

    dir_list.Delete();
}

#endif
