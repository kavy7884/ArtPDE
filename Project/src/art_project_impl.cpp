
ArtProjectBuilder ArtProject::create(const std::string & projectName) {
	return ArtProjectBuilder(projectName);
}
#if defined( _MSC_VER )

int ArtProject::checkFolder(const std::string &folderPath) {
	char sep = *this->slash.begin();

	//std::cout << sep << std::endl;

	for (std::string::const_iterator iter = folderPath.cbegin(); iter != folderPath.cend(); )
	{
		std::string::const_iterator newIter = std::find(iter, folderPath.cend(), sep);

		if (newIter == folderPath.cbegin()) {
			++iter;
			newIter = std::find(iter, folderPath.cend(), sep);
		}

		std::string newPath = std::string(folderPath.cbegin(), newIter);

		DWORD dwAttrib = GetFileAttributesA(newPath.c_str());

		if (!(dwAttrib != INVALID_FILE_ATTRIBUTES &&
			(dwAttrib & FILE_ATTRIBUTE_DIRECTORY)))
		{
			if (_mkdir(newPath.c_str()) != 0)
			{
				std::cout << ">> cannot create folder [" << newPath << "] ! " << std::endl;
				return -1;
			}
			else
				std::cout << ">> folder [" << newPath << "] created! " << std::endl;
		}
	
		iter = newIter;
		if (newIter != folderPath.end())
			++iter;
	}
	return 0;
}

#else
int ArtProject::checkFolder(const std::string &folderPath) {

	struct stat st;
	mode_t mode = 0777;
	char sep = *this->slash.begin();

	for (std::string::const_iterator iter = folderPath.cbegin(); iter != folderPath.cend(); )
	{
		std::string::const_iterator newIter = std::find(iter, folderPath.cend(), sep);

		if (newIter == folderPath.cbegin()) {
			++iter;
			newIter = std::find(iter, folderPath.cend(), sep);
		}

		std::string newPath = std::string(folderPath.cbegin(), newIter);

		if (stat(newPath.c_str(), &st) != 0)
		{
			if (mkdir(newPath.c_str(), mode) != 0 && errno != EEXIST)
			{
				std::cout << ">> cannot create folder [" << newPath << "] : " << strerror(errno) << std::endl;
				return -1;
			}
			else
				std::cout << ">> folder [" << newPath << "] created! " << std::endl;

		}
		else
			if (!S_ISDIR(st.st_mode))
			{
				errno = ENOTDIR;
				std::cout << ">> path [" << newPath << "] not a dir " << std::endl;
				return -1;
			}
		//        else
		//            std::cout << ">> path [" << newPath << "] already exists " << std::endl;

		iter = newIter;
		if (newIter != folderPath.end())
			++iter;
	}
	return 0;
}
#endif

bool ArtProject::Init() {
	if (checkFolder(getProjectPath()) != 0) return false;
	if (checkFolder(getProjectGeometryPath()) != 0) return false;
	if (checkFolder(getProjectInitialPath()) != 0) return false;
	if (checkFolder(getProjectSettingPath()) != 0) return false;
	if (checkFolder(getProjectResultsPath()) != 0) return false;

	return true;
}